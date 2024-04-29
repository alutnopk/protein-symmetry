#!/bin/bash
# - Biological Assembly File in PDB: `wget https://files.rcsb.org/download/1hh3.pdb1`
# - Biological Assembly File in PDBx/mmCIF: `wget https://files.rcsb.org/download/5a9z-assembly1.cif`
# - PDB: `wget https://files.rcsb.org/download/4hhb.pdb`
# - PDBx/mmCIF: `wget https://files.rcsb.org/download/4hhb.cif`

if [ -z "$1" ]; then
    echo "Enter an argument."
    exit 1
fi

wget https://files.rcsb.org/download/"$1".gz && gunzip "$1".gz