import pymol
from pymol import cmd
from ProteinSymmetry import *
import numpy as np
import sys

AXIS_LEN = 75

def main():
	# TODO: render axes using CGO objects
	cmd.remove("hetatm")
	cmd.util.cbc()
	C, subC = mark_centroids()
	join_centroids(subC)
	n = len(subC)
	print(n)
	# find partitions
	partition1 = {}
	partition2 = {}
	partition1, partition2 = create_partition(subC)
	# recurse(subC.keys(), is_regular_polygon, subC)

	# primary axis
	prim1 = np.mean([subC[chain] for chain in partition1], axis=0)
	prim2 = np.mean([subC[chain] for chain in partition2], axis=0)
	prim_axis = (prim1 - prim2) / np.linalg.norm(prim1 - prim2)
	cmd.pseudoatom("prim1", pos=list(C + AXIS_LEN * prim_axis))
	cmd.pseudoatom("prim2", pos=list(C - AXIS_LEN * prim_axis))
	cmd.distance("prim_axis", "prim1", "prim2")
	# secondary axis
	chain1 = "A"
	max_dist = 0.0
	nbr = ''
	for chain2 in partition2:
		dist = np.linalg.norm(subC[chain1] - subC[chain2])
		if dist > max_dist:
			max_dist = dist
			nbr = chain2
	sec_axis = np.cross(subC[chain1] - subC[nbr], prim_axis)
	sec_axis /= np.linalg.norm(sec_axis)
	cmd.pseudoatom(f"sec11", pos=list(C + AXIS_LEN * sec_axis))
	cmd.pseudoatom(f"sec12", pos=list(C - AXIS_LEN * sec_axis))
	cmd.distance(f"sec_axis_1", f"sec11", f"sec12")
	# remaining secondary axes can be found easily via rotation
	# rotate sec_axis about prim_axis by pi / n

	# using euler rodrigues formula
	theta = 2 * np.pi / n
	a = np.cos(theta / 2.0)
	b, c, d = - prim_axis * np.sin(theta / 2.0)
	rot_mat = np.array([[a*a + b*b - c*c - d*d, 2 * (b*c + a*d), 2 * (b*d - a*c)],
					 [2 * (b*c - a*d), a*a + c*c - b*b - d*d, 2 * (c*d + a*b)],
					 [2 * (b*d + a*c), 2 * (c*d - a*b), a*a + d*d - b*b - c*c]])
	sec_axis_new = np.copy(sec_axis)
	# iterate from pi/n to (n-1)pi/n
	for i in range(1, n):
		sec_axis_new = rot_mat @ sec_axis_new
		print(sec_axis_new)
		sec_axis_new /= np.linalg.norm(sec_axis_new)
		cmd.pseudoatom(f"sec{i+1}1", pos=list(C + AXIS_LEN * sec_axis_new))
		cmd.pseudoatom(f"sec{i+1}2", pos=list(C - AXIS_LEN * sec_axis_new))
		cmd.distance(f"sec_axis_{i+1}", f"sec{i+1}1", f"sec{i+1}2")

	# # quaternion approach
	# angl = 2 * np.pi / n
	# Q = np.array([	np.cos(angl/2), 
	# 				np.sin(angl/2) * prim_axis[0],
	# 				np.sin(angl/2) * prim_axis[1],
	# 				np.sin(angl/2) * prim_axis[2]
	# 			])
	# Q /= np.linalg.norm(Q)
	# Q_inv = np.array([Q[0], -Q[1], -Q[2], -Q[3]])
	# P = np.array([0.0, sec_axis[0], sec_axis[1], sec_axis[2]])
	# R = (Q @ P) @ Q_inv
	# sec_axis_new = R[1:]

if __name__ == '__main__':
	main()
