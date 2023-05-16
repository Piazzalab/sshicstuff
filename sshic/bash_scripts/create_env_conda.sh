#!/bin/bash

# Check if the yml file exists
if [ $# -eq 0 ]; then
    echo "Veuillez fournir le chemin vers le fichier YAML."
    exit 1
fi

# Path to yml file
yaml_file="$1"

# conda env's name
env_name=$(basename "$yaml_file" .yml)

# creation of the env
conda env create --name "$env_name" --file "$yaml_file"

# Vcheck if the creation has succeed
if [ $? -eq 0 ]; then
    echo "L'environnement Conda '$env_name' a été créé avec succès."
else
    echo "Une erreur s'est produite lors de la création de l'environnement Conda."
fi

