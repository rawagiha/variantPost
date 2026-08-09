FROM python:3.10-slim

RUN apt-get update && apt-get install -y \
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

RUN pip install --no-cache-dir Cython pysam numpy pandas scipy

COPY . /app
RUN pip install --no-cache-dir .

ENTRYPOINT ["indelinside"]
