from Bio.PDB import *
import pymol
from pymol import cmd

def cyclic_axes(C_main, C_sub):
	'''Calculate candidate axis for cyclic symmetry as a numpy array, given the centroid and subcentroids
	'''
	n = len(C_sub)
	# C_1, C_2 = C_sub[0], C_sub[1]
	# cyclic_axis = np.cross(C_1 - C_main, C_2 - C_main)
	cyclic_axis = np.zeros(3)
	for i in range(n):
		for j in range(i+1, n):
			C_i, C_j = C_sub[i], C_sub[j]
			cyclic_axis += np.cross(C_i - C_main, C_j - C_main)
	cyclic_axis /= np.linalg.norm(cyclic_axis, 2)
	return cyclic_axis

def fetch_axes(entry_id, assembly_id=1):
	'''Fetch the annotated axes from RCSB PDB'''
	axes = []
	import requests as rq, json, sys
	resp = rq.get(f"https://data.rcsb.org/rest/v1/core/assembly/{entry_id}/{assembly_id}")
	# resp = rq.get(f"https://data.rcsb.org/rest/v1/core/polymer_entity/{entry_id}/{assembly_id}")
	if resp.status_code == 200:
		data = resp.json()
		if 'rcsb_struct_symmetry' in data:
			print(json.dumps(data, indent=4))
			# print(json.dumps(data['rcsb_struct_symmetry'], indent=4))
			for rot_axis in data['rcsb_struct_symmetry'][0]['rotation_axes']:
				axes.append([rot_axis['start'], rot_axis['end']])
		else:
			print("'rcsb_struct_symmetry' not found")
	else:
		print(f"Response {resp.status_code}", file=sys.stderr)
	# print(axes)
	print("------------------------------------------------------------------")
	resp = rq.get(f"https://data.rcsb.org/rest/v1/core/entry/{entry_id}")
	data = resp.json()
	print(json.dumps(data, indent=4))
	return axes

def load_protein(pdb_file):
	'''Start a Pymol process and load up a protein via PDB ID'''
	pymol.finish_launching(["pymol"])
	# cmd.fetch(pdb_id, name=pdb_id, type='pdb1')
	cmd.load(pdb_file)

def display_axes(axes):
	'''Render the desired axes'''
	for i in range(0, len(axes)):
		axis = axes[i]
		cmd.pseudoatom(f'A{i}1', pos=axis[0], color='yellow')
		cmd.pseudoatom(f'A{i}2', pos=axis[1], color='yellow')
		cmd.distance(f'axis{i}', f'A{i}1', f'A{i}2')
	cmd.orient()

def print_structure(structure):
	# print(f"structure centred at: {structure.center_of_mass()}")
	# for model in structure:
	# 	print(f"\tmodel centred at: {model.center_of_mass()}")
	# 	for chain in model:
	# 		print(f"\t\tchain centred at: {chain.center_of_mass()}")
	# 		for residue in chain:
	# 			# print(f"\t\t\tresidue centred at: {residue.center_of_mass()}")
	# 			for atom in residue:
	# 				# print(atom)
	pass