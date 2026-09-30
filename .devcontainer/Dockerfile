# Bioinformatics teaching image for GitHub Codespaces
# Ubuntu + Miniforge (conda/mamba) + one "bioinfo" environment with the whole toolkit.
FROM mcr.microsoft.com/devcontainers/base:ubuntu-24.04

ARG USERNAME=vscode
ENV CONDA_DIR=/opt/conda \
    ENV_NAME=bioinfo
# Put the bioinfo env first on PATH so tools work in every shell, including
# the non-interactive shells Nextflow uses to run tasks.
ENV PATH=${CONDA_DIR}/envs/${ENV_NAME}/bin:${CONDA_DIR}/bin:${PATH}
# Keep the Nextflow JVM modest on an 8 GB machine
ENV NXF_OPTS="-Xms256m -Xmx2g"

RUN apt-get update \
 && export DEBIAN_FRONTEND=noninteractive \
 && apt-get install -y --no-install-recommends wget ca-certificates bzip2 procps tree less \
 && apt-get clean && rm -rf /var/lib/apt/lists/*

# Install conda as the vscode user so students can add packages without sudo
RUN mkdir -p ${CONDA_DIR} && chown ${USERNAME}:${USERNAME} ${CONDA_DIR}
USER ${USERNAME}

RUN wget -q "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-$(uname -m).sh" -O /tmp/miniforge.sh \
 && bash /tmp/miniforge.sh -b -u -p ${CONDA_DIR} \
 && rm /tmp/miniforge.sh

COPY --chown=${USERNAME}:${USERNAME} environment.yml /tmp/environment.yml
RUN mamba env create -f /tmp/environment.yml \
 && mamba clean -afy \
 && rm /tmp/environment.yml

# Channel setup so "mamba install <tool>" just works, auto-activate the env,
# and warm up Nextflow so the first run in class is fast.
RUN conda config --add channels bioconda \
 && conda config --add channels conda-forge \
 && conda config --set channel_priority strict \
 && conda init bash \
 && echo "conda activate ${ENV_NAME}" >> ~/.bashrc \
 && (nextflow info > /dev/null 2>&1 || true)

USER root
