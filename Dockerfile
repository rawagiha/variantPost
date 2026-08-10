FROM mambaorg/micromamba:1.5.8-git-apt

RUN micromamba install -y -n base -c conda-forge -c bioconda \
    python=3.10 \
    cython \
    numpy \
    scipy \
    pysam \
    pandas \
    c-compiler \
    cxx-compiler \
    && micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1

WORKDIR /app
COPY . /app

RUN pip install --no-cache-dir --no-build-isolation .

ENTRYPOINT ["indelinside"]
