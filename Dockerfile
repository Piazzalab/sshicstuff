# =========================
# Stage 1 — build oligo4sshic
# =========================
FROM rust:1.81-bullseye AS o4s-builder

# If/when oligo4sshic is public, update the URL below (HTTPS clone).
# You can also override at build time:
#   docker build --build-arg O4S_GIT_URL=https://github.com/.../oligo4sshic.git -t sshicstuff:latest .
ARG O4S_GIT_URL=https://gitbio.ens-lyon.fr/LBMC/GM/oligo4sshic.git

WORKDIR /src
RUN git clone --depth=1 "$O4S_GIT_URL" oligo4sshic \
 && cd oligo4sshic \
 && cargo build --release

# The binary ends up at: /src/oligo4sshic/target/release/oligo4sshic


# =========================
# Stage 2 — app + conda env
# =========================
FROM mambaorg/micromamba:1.5.10

# ---- System packages (lean) ----
USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
      build-essential git ca-certificates curl \
      zlib1g-dev libbz2-dev liblzma-dev libssl-dev libffi-dev \
  && rm -rf /var/lib/apt/lists/*

# ---- Workdir & copy metadata first (cache-friendly) ----
WORKDIR /app
COPY environment.yml .
COPY requirements.txt .
COPY pyproject.toml .
COPY README.md .

# ---- Create the conda env from environment.yml ----
# This typically brings in CLI tools (bowtie2, minimap2, samtools, seqtk, pairtools, etc.).
ENV MAMBA_DOCKERFILE_ACTIVATE=1
SHELL ["/bin/bash", "-lc"]
RUN micromamba create -y -n sshicstuff_env -f environment.yml \
 && micromamba clean -a -y

# ---- Python deps (pin/upgrade as needed) ----
RUN micromamba run -n sshicstuff_env python -m pip install --upgrade pip \
 && micromamba run -n sshicstuff_env python -m pip install -r requirements.txt

# ---- App source & install (editable for dev; switch to non-editable for release) ----
COPY . /app
RUN micromamba run -n sshicstuff_env python -m pip install -e .

# ---- Bring in the oligo4sshic binary from the builder stage ----
COPY --from=o4s-builder /src/oligo4sshic/target/release/oligo4sshic /usr/local/bin/oligo4sshic

# ---- Non-root user & PATH ----
RUN useradd -ms /bin/bash appuser && chown -R appuser:appuser /app
USER appuser
ENV PATH=/opt/conda/envs/sshicstuff_env/bin:$PATH

# If you run the Dash/Flask UI (e.g., `sshicstuff view`), 8050 is common.
EXPOSE 3838

# The CLI name exposed by your pyproject's [project.scripts] -> "sshicstuff = sshicstuff.__main__:main"
ENTRYPOINT ["sshicstuff"]
CMD ["view"]
