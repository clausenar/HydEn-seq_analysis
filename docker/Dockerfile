FROM continuumio/miniconda3:latest

# Pin to a branch, tag, or full commit SHA for reproducible builds.
ARG GIT_REF=main
ARG REPO_URL=https://github.com/clausenar/HydEn-seq_analysis.git

# git + AWS CLI (git clones the pipeline itself below; awscli drives the S3
# sync in batch_entrypoint.sh)
RUN apt-get update && apt-get install -y --no-install-recommends curl unzip git \
    && curl -sL "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o /tmp/awscliv2.zip \
    && unzip -q /tmp/awscliv2.zip -d /tmp \
    && /tmp/aws/install \
    && rm -rf /tmp/awscliv2.zip /tmp/aws \
    && apt-get purge -y unzip && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*

# Avoid Anaconda's "defaults" channel entirely (older conda in this base
# image doesn't support the "conda tos accept" ToS flow it now requires) --
# conda-forge/bioconda are used for everything below and need no ToS.
RUN conda config --remove channels defaults || true \
    && conda config --add channels conda-forge

# bio_env: analysis/plotting (count_bases.py, origin_metaplot.py)
RUN conda create -n bio_env python=3.10 -y \
    && conda install -n bio_env -c bioconda -c conda-forge pybedtools pysam bedtools samtools -y \
    && conda install -n bio_env -c conda-forge matplotlib -y \
    && conda clean -afy

# hyden_pipeline: raw read processing (runs the Snakefile itself)
RUN conda create -n hyden_pipeline -c bioconda -c conda-forge python=3.11 snakemake bowtie samtools bedtools perl-dbi perl-dbd-mysql seqkit -y \
    && conda run -n hyden_pipeline pip install cutadapt \
    && conda clean -afy

# Clone the pipeline itself - the image runs the exact, unmodified
# git-tracked code at GIT_REF, not a hand-patched local copy. Works for a
# branch, tag, or full commit SHA.
RUN git init /app \
    && cd /app \
    && git remote add origin "${REPO_URL}" \
    && git fetch --depth 1 origin "${GIT_REF}" \
    && git checkout FETCH_HEAD \
    && chmod +x /app/docker/batch_entrypoint.sh

WORKDIR /app
ENTRYPOINT ["/app/docker/batch_entrypoint.sh"]
