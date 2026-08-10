FROM python:3.10-slim

RUN rm -f /etc/system-fips /etc/crypto-policies/backends/openssl.config

ENV OPENSSL_FORCE_FIPS_MODE=0
ENV OPENSSL_FIPS=0
ENV OPENSSL_CONF=/dev/null

RUN apt-get update && apt-get install -y \
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

RUN pip install --no-cache-dir setuptools wheel Cython pysam numpy pandas scipy

COPY . /app

RUN pip install --no-cache-dir --no-build-isolation .

ENTRYPOINT ["indelinside"]
