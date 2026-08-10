FROM python:3.10-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    git \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

RUN pip install --no-cache-dir --upgrade pip setuptools wheel
RUN pip install --no-cache-dir Cython numpy scipy pysam pandas

COPY . /app
RUN pip install --no-cache-dir --no-build-isolation .

ENV OPENSSL_FORCE_FIPS_ZERO=1
ENV OPENSSL_CONF=/dev/null

ENTRYPOINT ["indelinside"]
