#!/bin/bash
# Argument types:
# - Biological Assembly File in PDB: 1hh3.pdb1
# - Biological Assembly File in PDBx/mmCIF: 5a9z-assembly1.cif
# - PDB: 4hhb.pdb
# - PDBx/mmCIF: 1a0d.cif

if [ -z "$1" ]; then
    echo "Enter the full name of the resource, with the extension."
    exit 1
fi

wget https://files.rcsb.org/download/"$1".gz && gunzip "$1".gz