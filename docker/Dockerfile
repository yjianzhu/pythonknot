FROM quay.io/pypa/manylinux2014_x86_64:latest

ENV RUSTUP_HOME=/opt/rust/rustup
ENV CARGO_HOME=/opt/rust/cargo
ENV PATH=/opt/rust/cargo/bin:${PATH}

RUN curl https://sh.rustup.rs -sSf | sh -s -- -y --profile minimal \
    && rustc --version \
    && cargo --version \
    && yum install -y patchelf \
    && /opt/python/cp313-cp313/bin/python -m pip install -U auditwheel \
    && /opt/python/cp313-cp313/bin/python -m auditwheel --version \
    && yum clean all \
    && rm -rf /var/cache/yum

WORKDIR /work/pythonknot
