#!/bin/bash
# Setup script for ancestral sequence reconstruction pipeline
# Usage: bash setup.sh

set -e

# Create conda environment if it doesn't exist
env_name="ancestral_reconstruction"
if ! conda info --envs | grep -q "^$env_name"; then
    echo "Creating conda environment $env_name..."
    conda env create -f env/env.yaml
else
    echo "Conda environment $env_name already exists."
fi

echo "Activating environment and running pipeline..."
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate $env_name
nextflow run main.nf
