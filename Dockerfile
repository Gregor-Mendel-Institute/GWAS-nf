FROM mambaorg/micromamba:0.9.2
LABEL authors="p.huether@lmu.de" \
    description="Container image containing all dependencies for the GWAS-nf pipeline"

COPY environment.yml /environment.yml

RUN micromamba install -y -n base -f /environment.yml && micromamba clean -a
