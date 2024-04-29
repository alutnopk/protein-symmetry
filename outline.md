# Goal
Figure out a complete algorithm for D2, D3, D4

# Tasks

- [] Learn pymol object selection
- [] Try using cmd.align and related commands
- [] Management scripts:
	- fetch and render the annotated axes from RCSB
	- quickly drawing axes and showing subcentroids
	- quick verification of accuracy
- [] Start by looking at C1, C2, C4 and D2 homo-tetramers and identify how to separate them
- [] Move on to heteromers
- [] Collect dataset and test the algorithm

- Biological Assembly File in PDB: `wget https://files.rcsb.org/download/1hh3.pdb1`
- Biological Assembly File in PDBx/mmCIF: `wget https://files.rcsb.org/download/5a9z-assembly1.cif`
- PDB: `wget https://files.rcsb.org/download/4hhb.pdb`
- PDBx/mmCIF: `wget https://files.rcsb.org/download/4hhb.cif`
- `rcsbsearchapi`

## Pymol commands
- pseudoatom my_atom, pos=[10.0, 10.0, 10.0]
- distance cyclic_axis, atom1, atom2
- remember to extend the line

## Biopython source files
- PDBParser
- StructureBuilder
- Structure
- Entity
- vectors
- Chain