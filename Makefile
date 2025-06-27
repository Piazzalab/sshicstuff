# Variables
ENV_NAME := sshicstuff_env
YAML := environment.yml
PYTHON := $(shell which python)

.PHONY: all env install-oligo4sshic clean

all: env install install-oligo4sshic

# Step 1: Create conda environment
env:
	@echo "Creating Conda environment from $(YAML)..."
	@mamba env create -f $(YAML) -n $(ENV_NAME)
	@echo "Environment $(ENV_NAME) created/updated."

# Step 2: pip install -e . in conda env
install:
	@echo "Installing sshicstuff in editable mode..."
	@mamba run -n $(ENV_NAME) pip install -e .
	@echo "Installed with pip editable mode."

# Step 3: install oligo4sshic
install-oligo4sshic:
	@echo "Checking if cargo is installed..."
	@if ! command -v cargo > /dev/null; then \
		echo "Rust not found, installing via rustup..."; \
		curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y; \
		source $$HOME/.cargo/env; \
	fi
	@echo "Cloning and installing oligo4sshic..."
	@git clone git@gitbio.ens-lyon.fr:LBMC/GM/oligo4sshic.git
	@cargo install --path oligo4sshic
	@echo "oligo4sshic installed."

# Cleanup
clean:
	@echo "Removing conda environment $(ENV_NAME)..."
	@conda remove -n $(ENV_NAME) --all -y
	@rm -rf oligo4sshic

