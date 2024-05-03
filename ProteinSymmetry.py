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
	pass

def show_pdb_axes(entry_id, assembly_id=1):
	# fetch the axes
	import requests as rq, json, sys
	pdb_axes = []
	resp = rq.get(f"https://data.rcsb.org/rest/v1/core/assembly/{entry_id}/{assembly_id}")
	if resp.status_code == 200:
		data = resp.json()
		if 'rcsb_struct_symmetry' in data:
			# print(json.dumps(data, indent=4))
			for rot_axis in data['rcsb_struct_symmetry'][0]['rotation_axes']:
				pdb_axes.append([rot_axis['start'], rot_axis['end']])
		else:
			print("'rcsb_struct_symmetry' not found")
	else:
		print(f"Response {resp.status_code}", file=sys.stderr)
	# # Debug printing
	# print("------------------------------------------------------------------")
	# resp = rq.get(f"https://data.rcsb.org/rest/v1/core/entry/{entry_id}")
	# data = resp.json()
	# print(json.dumps(data, indent=4))
	
	# display the axes
	for i in range(0, len(pdb_axes)):
		axis = pdb_axes[i]
		cmd.pseudoatom(f'P{i}1', pos=axis[0], color='white')
		cmd.pseudoatom(f'P{i}2', pos=axis[1], color='white')
		cmd.distance(f'pdb_axis{i}', f'P{i}1', f'P{i}2')

def mark_centroids(selection="all", center=0, quiet=1):
	chains = cmd.get_chains(selection)
	subcentr = {}

	for chain in chains:
		model = cmd.get_model(f"chain {chain}")
		print(type(model))
		centr = np.zeros(3)
		for atm in model.atom:
			centr += atm.coord
		centr /= len(model.atom)
		print(f"Centroid of chain {chain}: {centr}")
		cmd.pseudoatom("centr"+chain, pos=list(centr))
		subcentr[chain] = centr
	
	model = cmd.get_model(selection)
	centr = np.zeros(3)
	for atm in model.atom:
		centr += atm.coord
	centr /= len(model.atom)
	print(f"Centroid of assembly: {centr}")
	cmd.pseudoatom("centroid", pos=list(centr))
	# cmd.extend("mk", mark_centroids)
	return centr, subcentr

def join_centroids(subcentroids):
	for chain1 in subcentroids.keys():
		for chain2 in subcentroids.keys():
			if chain1 == chain2:
				continue
			cmd.distance(f"d{chain1}{chain2}", f"centr{chain1}", f"centr{chain2}")

def is_regular_polygon(chains, subC):
	# must be in same plane
	# must have same distance from centroid
	# each distance vector must be equally inclined
	# TODO
	centr = np.mean([subC[chain] for chain in chains], axis=0)
	np.sum([subC[chain] for chain in chains], axis=0)
	vec_sum = np.zeros(3)
	for chain in chains:
		vec_sum += centr - subC[chain]
	
	print(f"Error: {err}")
	if err <= 0.1:
		return True
	else:
		return False

def recurse(fn, subC, a = [], choice = [], i = 0):
	if a == []:
		a = subC.keys()
	if i == len(a):
		if len(choice) == len(a) / 2:
			fn(choice, subC)
		return
	choice.append(a[i])
	recurse(fn, subC, a, choice, i+1)
	choice.pop()
	recurse(fn, subC, a, choice, i+1)

def fetch_protein(pdb_id):
	import subprocess
	cmd_string = f"bash utils/getprot.sh {pdb_id}-assembly1.cif"
	subprocess.Popen(cmd_string, stdout=subprocess.PIPE, shell=True)
	output, error = process.communicate()
	if error:
		print(f"Error: {error}")
	else:
		print(f"Output: {output.decode('utf-8')}")

def load_protein(pdb_file):
	'''Start a PyMOL process and load up a protein via PDB ID'''
	pymol.finish_launching(["pymol", "-p"])
	cmd.load(pdb_file)

def create_partition(subC):
	n = len(subC)
	if n == 4:
		return {'A', 'C'}, {'B', 'D'}
	elif n == 6:
		return {'A', 'C', 'E'}, {'B', 'D', 'F'}
	elif n == 8:
		return {'A', 'A-2', 'A-3', 'A-4'}, {'B', 'B-2', 'B-3', 'B-4'}