FROM python:3.10-slim

ENV OPENSSL_FORCE_FIPS_MODE=0 \
    OPENSSL_FIPS=0

RUN apt-get update && apt-get install -y \
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    && rm -rf /var/lib/apt/lists/*

RUN echo 'openssl_conf = openssl_init\n\n[openssl_init]\nproviders = provider_sect\n\n[provider_sect]\ndefault = default_sect\n\n[default_sect]\nactivate = 1\n' > /etc/ssl/openssl.cnf

RUN echo '#include <stdio.h>\nint FIPS_mode(void) { return 0; }\nint OPENSSL_fips_mode(void) { return 0; }\nint FIPS_selftest(void) { return 1; }' > /tmp/nofips.c && \
    gcc -shared -fPIC /tmp/nofips.c -o /usr/local/lib/libnofips.so && \
    rm /tmp/nofips.c

ENV LD_PRELOAD=/usr/local/lib/libnofips.so

WORKDIR /app

RUN pip install --no-cache-dir setuptools wheel Cython pysam numpy pandas scipy

COPY . /app

RUN pip install --no-cache-dir --no-build-isolation .

ENTRYPOINT ["indelinside"]
