import sys
import os


AMINOACIDS = ["VAL", "ASN", "GLY", "LEU", "ILE",
              "SER", "ASP", "LYS", "MET", "GLN",
              "TRP", "ARG", "ALA", "THR", "PRO",
              "PHE", "GLU", "HIS", "HIP", "TYR",
              "CYS"]

TER_CONSTR = 2.5

BACK_CONSTR = 0.5


class ConstraintBuilder(object):

	def __init__(self, pdb):
		self.pdb = pdb

	def parse_atoms(self):

		protein_residues = 0
                initial_res = False

		with open(self.pdb, "r") as pdb:
			for line in pdb:
				if line.startswith("ATOM") and line[16:21].strip(' ') in AMINOACIDS:
					if not initial_res : initial_res = line[22:26].strip()
					else: protein_residues = line[22:26].strip()

		return int(initial_res), int(protein_residues)

	def build_constraint(self, initial_res, num_residues):

		constraints = ['''"constraints":[''', ]

		constraint_str = '''{{ "type": "constrainAtomToPosition", "springConstant": {1}, "equilibriumDistance": 0.0, "constrainThisAtom": "A:{0}:_CA_" }},'''

		for i in range(initial_res, num_residues, 10):
			if(i==1):
				constraint = constraint_str.format(i, TER_CONSTR)
			else:
				constraint = constraint_str.format(i, BACK_CONSTR)

			constraints.append(constraint)

		constraints.append(constraint_str.format(num_residues, TER_CONSTR).strip(","))

		constraints.append("],")

		return constraints

def retrieve_constraints(pdb_file):
	constr = ConstraintBuilder(pdb_file)
	initial_res, final_res = constr.parse_atoms()
        print(initial_res, final_res)
	constraints = constr.build_constraint(initial_res, final_res)
	return constraints




