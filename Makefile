SHELL := /bin/bash
.SHELLFLAGS := -eu -o pipefail -c

ENV_NAME    := sshicstuff_env
YAML        := environment.yml
LOCK_FILE   := conda-lock.yml
OLIGO_DIR   := oligo4sshic

MAMBA       ?= mamba
CONDA_LOCK  ?= conda-lock
CARGO       ?= cargo
GIT         ?= git

.PHONY: help all env env-update env-lock lock install-oligo4sshic reinstall clean clean-oligo clean-lock

help:
	@echo "Available targets:"
	@echo "  make env               Create the conda environment from $(YAML)"
	@echo "  make env-update        Update the existing conda environment from $(YAML)"
	@echo "  make lock              Generate $(LOCK_FILE) for linux-64 and osx-arm64"
	@echo "  make env-lock          Create the conda environment from $(LOCK_FILE)"
	@echo "  make reinstall         Reinstall the package in editable mode inside the env"
	@echo "  make install-oligo4sshic  Clone and install the Rust submodule"
	@echo "  make all               Create env + install oligo4sshic"
	@echo "  make clean             Remove the conda environment"
	@echo "  make clean-oligo       Remove cloned oligo4sshic directory"
	@echo "  make clean-lock        Remove $(LOCK_FILE)"

all: env install-oligo4sshic

# Create the conda environment from environment.yml.
# The editable install is triggered by the pip block in environment.yml.
env:
	@command -v $(MAMBA) >/dev/null 2>&1 || { echo "$(MAMBA) not found. Please install mamba first."; exit 1; }
	@if $(MAMBA) env list | awk '{print $$1}' | grep -qx "$(ENV_NAME)"; then \
		echo "Environment '$(ENV_NAME)' already exists. Use 'make env-update' or 'make clean' first."; \
		exit 1; \
	fi
	$(MAMBA) env create -f $(YAML)

# Update an existing environment from environment.yml.
env-update:
	@command -v $(MAMBA) >/dev/null 2>&1 || { echo "$(MAMBA) not found. Please install mamba first."; exit 1; }
	$(MAMBA) env update -n $(ENV_NAME) -f $(YAML) --prune

# Generate a fully pinned lock file for supported platforms.
lock:
	@command -v $(CONDA_LOCK) >/dev/null 2>&1 || { echo "$(CONDA_LOCK) not found. Run: pip install conda-lock"; exit 1; }
	$(CONDA_LOCK) lock \
		--file $(YAML) \
		--platform linux-64 \
		--lockfile $(LOCK_FILE)

# Create the environment from the lock file.
# This is the recommended fully reproducible installation path.
env-lock:
	@command -v $(CONDA_LOCK) >/dev/null 2>&1 || { echo "$(CONDA_LOCK) not found. Run: pip install conda-lock"; exit 1; }
	@test -f $(LOCK_FILE) || { echo "$(LOCK_FILE) not found. Run 'make lock' first."; exit 1; }
	@if $(MAMBA) env list | awk '{print $$1}' | grep -qx "$(ENV_NAME)"; then \
		echo "Environment '$(ENV_NAME)' already exists. Use 'make clean' first if you want a fresh locked install."; \
		exit 1; \
	fi
	$(CONDA_LOCK) install -n $(ENV_NAME) $(LOCK_FILE)

# Optional: force a local editable reinstall inside the environment.
# Useful during development, but not required for locked installs.
reinstall:
	@command -v $(MAMBA) >/dev/null 2>&1 || { echo "$(MAMBA) not found. Please install mamba first."; exit 1; }
	mamba run -n $(ENV_NAME) pip install -e .

# Clone and install the Rust submodule/binary.
install-oligo4sshic:
	@command -v $(GIT) >/dev/null 2>&1 || { echo "$(GIT) not found. Please install git first."; exit 1; }
	@if ! command -v $(CARGO) >/dev/null 2>&1; then \
		echo "cargo not found."; \
		echo "Please install Rust first: https://rustup.rs"; \
		exit 1; \
	fi
	@test -d $(OLIGO_DIR) || $(GIT) clone git@gitbio.ens-lyon.fr:LBMC/GM/oligo4sshic.git $(OLIGO_DIR)
	$(CARGO) install --path $(OLIGO_DIR)

clean:
	-$(MAMBA) env remove -n $(ENV_NAME) -y

clean-oligo:
	rm -rf $(OLIGO_DIR)

clean-lock:
	rm -f $(LOCK_FILE)