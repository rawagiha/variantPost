FROM condaforge/miniforge3:latest

RUN mamba install -y -c conda-forge -c bioconda \
    git \
    python=3.10 \
    cython \
    numpy \
    scipy \
    pysam \
    pandas \
    c-compiler \
    cxx-compiler \
    && mamba clean --all --yes

WORKDIR /app
COPY . /app

RUN pip install --no-cache-dir --no-build-isolation .

ENTRYPOINT ["indelinside"]
