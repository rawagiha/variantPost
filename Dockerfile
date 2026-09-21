FROM condaforge/miniforge3:latest

ENV CONDA_PKGS_DIRS=/tmp/pkgs

RUN mamba create -n app -y -c conda-forge -c bioconda \
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

ENV PATH=/opt/conda/envs/app/bin:$PATH

WORKDIR /app
COPY . /app

RUN pip install --no-cache-dir --no-build-isolation .

ENTRYPOINT ["indelinside"]
