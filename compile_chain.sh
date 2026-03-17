#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${ROOT_DIR}"

if [[ ! -d vendor/rust_knot ]]; then
  echo "ERROR: vendor/rust_knot not found. Run 'git submodule update --init --recursive' on host first."
  exit 1
fi

if ! command -v cargo >/dev/null 2>&1; then
  echo "ERROR: cargo not found in container. Install Rust toolchain first."
  exit 1
fi

if [[ -n "${PY_TAGS:-}" ]]; then
  read -r -a PY_TAGS_ARR <<< "${PY_TAGS}"
else
  PY_TAGS_ARR=(
    cp38-cp38
    cp39-cp39
    cp310-cp310
    cp311-cp311
    cp312-cp312
    cp313-cp313
    cp314-cp314
  )
fi

mkdir -p dist
rm -f dist/*.whl || true

for tag in "${PY_TAGS_ARR[@]}"; do
  if [[ ! -x "/opt/python/${tag}/bin/python" ]]; then
    echo "Skip ${tag}: interpreter not found"
    continue
  fi

  PYBIN="/opt/python/${tag}/bin"
  echo "=== Building wheel for ${tag} ==="
  rm -rf build
  "${PYBIN}/python" -m ensurepip --upgrade >/dev/null 2>&1 || true
  "${PYBIN}/python" -m pip install -U pip setuptools wheel
  "${PYBIN}/python" setup.py bdist_wheel
done

echo "Done. Raw wheels are in dist/."
