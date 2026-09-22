FROM condaforge/miniforge3:latest

RUN mamba create -n app -y -c conda-forge -c bioconda \
        python=3.10 \
        "pysam>=0.23.3" \
        "cython>=3.0.0" \
        numpy \
        scipy \
        pandas \
        c-compiler \
        cxx-compiler \
    && mamba clean --all --f --yes \
    && rm -rf /opt/conda/pkgs/* /tmp/*

ENV PATH=/opt/conda/envs/app/bin:$PATH

WORKDIR /app

COPY . /app

RUN rm -rf build/ dist/ *.egg-info variantpost/*.so \
    && pip install --no-cache-dir --no-build-isolation .

ENTRYPOINT ["indelinside"]
