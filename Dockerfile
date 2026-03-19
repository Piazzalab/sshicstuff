# Stage 1 — build the oligo4sshic binary
FROM rust:1.81-bullseye AS o4s-builder

# Override at build time if the repository becomes public:
#   docker build --build-arg O4S_GIT_URL=https://... -t sshicstuff:latest .
ARG O4S_GIT_URL=https://gitbio.ens-lyon.fr/LBMC/GM/oligo4sshic.git

WORKDIR /src
RUN git clone --depth=1 "$O4S_GIT_URL" oligo4sshic \
 && cd oligo4sshic \
 && cargo build --release


# Stage 2 — conda environment and application
FROM mambaorg/micromamba:1.5.10

USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
      build-essential git ca-certificates curl \
      zlib1g-dev libbz2-dev liblzma-dev libssl-dev libffi-dev \
  && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Copy dependency files before the source so Docker can cache the env layer.
COPY environment.yml requirements.txt pyproject.toml README.md ./

# Create the conda environment.
# The pip block in environment.yml runs "pip install -e ." which reads
# requirements.txt through pyproject.toml, so no separate pip step is needed.
ENV MAMBA_DOCKERFILE_ACTIVATE=1
SHELL ["/bin/bash", "-lc"]
RUN micromamba create -y -n sshicstuff_env -f environment.yml \
    && micromamba clean -a -y

# Copy the full source after the env is built to preserve the cache layer above.
COPY . /app
RUN micromamba run -n sshicstuff_env pip install --no-deps -e .

COPY --from=o4s-builder /src/oligo4sshic/target/release/oligo4sshic /usr/local/bin/oligo4sshic

RUN useradd -ms /bin/bash appuser \
    && chown -R appuser:appuser /app \
    && mkdir -p /cache \
    && chown -R appuser:appuser /cache

USER appuser
ENV PATH=/opt/conda/envs/sshicstuff_env/bin:$PATH
ENV SSHICSTUFF_CACHE_DIR=/cache
ENV FLASK_SECRET_KEY="xOt9KYbBDN4Fm84bzq2tUhs9PXjN6tGH8j7s3R9zNaPpQWqs"

EXPOSE 8050
ENTRYPOINT ["sshicstuff"]
CMD ["view"]