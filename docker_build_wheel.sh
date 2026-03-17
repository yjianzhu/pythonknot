#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${ROOT_DIR}"

IMAGE_NAME="${IMAGE_NAME:-pythonknot_compile}"
CONTAINER_NAME="${CONTAINER_NAME:-pythonknot_build}"
PY_TAGS="${PY_TAGS:-cp314-cp314}"
REPAIR_WHEEL="${REPAIR_WHEEL:-1}"
REPAIR_PLAT="${REPAIR_PLAT:-manylinux2014_x86_64}"
AUDITWHEEL_PY_TAG="${AUDITWHEEL_PY_TAG:-cp313-cp313}"

docker build -t "${IMAGE_NAME}" -f docker/Dockerfile .

TARGET_IMAGE_ID="$(docker image inspect -f '{{.Id}}' "${IMAGE_NAME}")"
if ! docker container inspect "${CONTAINER_NAME}" >/dev/null 2>&1; then
  docker create \
    --name "${CONTAINER_NAME}" \
    -v "${ROOT_DIR}:/work/pythonknot" \
    "${IMAGE_NAME}" \
    tail -f /dev/null >/dev/null
else
  CONTAINER_IMAGE_ID="$(docker inspect -f '{{.Image}}' "${CONTAINER_NAME}")"
  if [[ "${CONTAINER_IMAGE_ID}" != "${TARGET_IMAGE_ID}" ]]; then
    docker rm -f "${CONTAINER_NAME}" >/dev/null
    docker create \
      --name "${CONTAINER_NAME}" \
      -v "${ROOT_DIR}:/work/pythonknot" \
      "${IMAGE_NAME}" \
      tail -f /dev/null >/dev/null
  fi
fi

if [[ "$(docker inspect -f '{{.State.Running}}' "${CONTAINER_NAME}")" != "true" ]]; then
  docker start "${CONTAINER_NAME}" >/dev/null
fi

docker exec \
  -e PY_TAGS="${PY_TAGS}" \
  -e REPAIR_WHEEL="${REPAIR_WHEEL}" \
  -e REPAIR_PLAT="${REPAIR_PLAT}" \
  -e AUDITWHEEL_PY_TAG="${AUDITWHEEL_PY_TAG}" \
  -e HOST_UID="$(id -u)" \
  -e HOST_GID="$(id -g)" \
  "${CONTAINER_NAME}" \
  /bin/bash -lc '
    set -euo pipefail
    cd /work/pythonknot
    ./compile_chain.sh
    if [[ "${REPAIR_WHEEL}" == "1" ]]; then
      if [[ ! -x "/opt/python/${AUDITWHEEL_PY_TAG}/bin/python" ]]; then
        echo "ERROR: /opt/python/${AUDITWHEEL_PY_TAG}/bin/python not found."
        exit 1
      fi
      if ! /opt/python/${AUDITWHEEL_PY_TAG}/bin/python -c "import auditwheel" >/dev/null 2>&1; then
        /opt/python/${AUDITWHEEL_PY_TAG}/bin/python -m pip install -U auditwheel patchelf
      fi
      export PATH="/opt/python/${AUDITWHEEL_PY_TAG}/bin:${PATH}"
      mkdir -p wheelhouse
      /opt/python/${AUDITWHEEL_PY_TAG}/bin/python -m auditwheel repair \
        --plat "${REPAIR_PLAT}" \
        -w wheelhouse \
        dist/*.whl
    fi
    chown -R "${HOST_UID}:${HOST_GID}" dist wheelhouse build src/pythonknot.egg-info || true
  '

if [[ "${REPAIR_WHEEL}" == "1" ]]; then
  echo "Wheel build complete."
  echo "Raw wheels: dist/"
  echo "Repaired wheels: wheelhouse/ (${REPAIR_PLAT})"
else
  echo "Wheel build complete. Raw wheels: dist/"
fi
echo "Container: ${CONTAINER_NAME}; tags: ${PY_TAGS}"
