FROM continuumio/miniconda3

RUN conda install -c conda-forge autodock-vina openbabel pandas -y && \
    conda clean -afy

WORKDIR /docking
