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

USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
      build-essential git ca-certificates curl \
      zlib1g-dev libbz2-dev liblzma-dev libssl-dev libffi-dev \
  && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# ---- Copy only dependency files (for caching) ----
COPY environment.yml requirements.txt pyproject.toml README.md /app/

# ---- Build frozen env (atomic, cachable) ----
ENV MAMBA_DOCKERFILE_ACTIVATE=1
SHELL ["/bin/bash", "-lc"]
RUN micromamba create -y -n sshicstuff_env -f environment.yml \
    && micromamba run -n sshicstuff_env pip install --upgrade pip \
    && micromamba run -n sshicstuff_env pip install -r requirements.txt \
    && micromamba clean -a -y

# ---- Copy the full app only AFTER ----
COPY . /app
RUN micromamba run -n sshicstuff_env pip install -e .

# === Reste inchangé ===
COPY --from=o4s-builder /src/oligo4sshic/target/release/oligo4sshic /usr/local/bin/oligo4sshic

RUN useradd -ms /bin/bash appuser && chown -R appuser:appuser /app
RUN mkdir -p /cache && chown -R appuser:appuser /cache

USER appuser
ENV PATH=/opt/conda/envs/sshicstuff_env/bin:$PATH
ENV SSHICSTUFF_CACHE_DIR=/cache
ENV FLASK_SECRET_KEY="xOt9KYbBDN4Fm84bzq2tUhs9PXjN6tGH8j7s3R9zNaPpQWqs"

EXPOSE 8050
ENTRYPOINT ["sshicstuff"]
CMD ["view"]