#!/bin/bash

# Name of the Conda environment (change to your environment name)
CONDA_ENV="AlphaPulldown2"

# Check if the input PDB file is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: ./convert_pdb.sh <input.pdb>"
    exit 1
fi

PDB_FILE=$1

# Activate Conda environment
source ~/.bashrc
conda activate $CONDA_ENV

# Run the Python script
python pdb_to_cif.py "$PDB_FILE"

# Deactivate Conda environment
conda deactivate