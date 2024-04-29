import pymol
from pymol import cmd
from ProteinSymmetry import *
import sys

pdb_file = sys.argv[1]
pdb_id = sys.argv[2]
load_protein(pdb_file)
rcsb_axes = fetch_axes(pdb_id)
print(rcsb_axes)
display_axes(rcsb_axes)