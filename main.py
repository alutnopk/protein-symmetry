from ProteinSymmetry import *
import pymol
from pymol import cmd
import sys

def main():
	if len(sys.argv) < 2:
		print("Usage: python3 <script> <pdb-file>", file=sys.stderr)
		sys.exit(1)
	# pymol.finish_launching(['pymol'])
	# cmd.fetch(sys.argv[1], sys.argv[1])
	parser = PDBParser() # or MMCIFParser()
	structure = parser.get_structure("FOO", sys.argv[1]) # TODO: add exception handling
	C_main = structure.center_of_mass(geometric=True)
	# structure -> model -> chain -> residue -> atom
	C_sub = [model.center_of_mass(geometric=True) for model in structure]

	if len(C_sub) >= 3:
		print(f"Center: {C_main}")
		print(f"Direction ratios: {cyclic_axes(C_main, C_sub)}")
		rcsb_axes = fetch_axes('1a0l')
		display_axes(rcsb_axes)
		# then transform structure about cyclic_axis and check RMSD

	elif len(C_sub) == 2:
		# Dimer
		print(f"Dimer case, will check orientation of dimers to verify symmetry")
	else:
		print("Single chain / Error parsing the PDB")

if __name__ == '__main__':
	main()