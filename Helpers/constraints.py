import sys
import os


AMINOACIDS = ["VAL", "ASN", "GLY", "LEU", "ILE",
              "SER", "ASP", "LYS", "MET", "GLN",
              "TRP", "ARG", "ALA", "THR", "PRO",
              "PHE", "GLU", "HIS", "HIP", "TYR",
              "CYS"]

TER_CONSTR = 2.5

BACK_CONSTR = 0.5

CONSTR_ATOM = '''{{ "type": "constrainAtomToPosition", "springConstant": {0}, "equilibriumDistance": 0.0, "constrainThisAtom": "{1}:{2}:{3}" }},'''

CONSTR_DIST = '''{{ "type": "constrainAtomsDistance", "springConstant": {}, "equilibriumDistance": {}, "constrainThisAtom": "{}:{}:{}", "toThisOtherAtom": "{}:{}:{}" }},'''

CONSTR_CALPHA = '''{{ "type": "constrainAtomToPosition", "springConstant": {1}, "equilibriumDistance": 0.0, "constrainThisAtom": "A:{0}:_CA_" }},'''

class ConstraintBuilder(object):

	def __init__(self, pdb, gaps, metals):
		self.pdb = pdb
                self.gaps = gaps
                self.metals = metals

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

                
		init_constr = ['''"constraints":[''', ]

		back_constr = [CONSTR_CALPHA.format(i, BACK_CONSTR) for i in range(initial_res+1, num_residues, 10)]

		gaps_constr = self.gaps_constraints()

		metal_constr = self.metal_constraints()		

                terminal_constr = [CONSTR_CALPHA.format(initial_res, TER_CONSTR), CONSTR_CALPHA.format(num_residues, TER_CONSTR).strip(",")]

                final_constr = [ "],"]

		constraints = init_constr + back_constr + gaps_constr + metal_constr + terminal_constr + final_constr

		return constraints
	
	def gaps_constraints(self):
                gaps_constr = []
       		for  chain, residues in self.gaps.iteritems():
                 	gaps_constr = [CONSTR_ATOM.format(10, chain, residue, "_CA_") for residue in residues]
		return gaps_constr
	
	def metal_constraints(self):
		metal_constr = []
		for metal, ligands in self.metals.iteritems():
			metal_name, chain, metnum = metal.split(" ")
                        for ligand in ligands:
				ligand_info, bond_lenght = ligand
				resnum, resname, chain, ligname = ligand_info.split(" ")
                        	metal_constr.append(CONSTR_DIST.format(5, bond_lenght, chain, resnum, ligname, chain, metnum, metal_name))				
		return metal_constr

        

def retrieve_constraints(pdb_file, gaps, metal):
	constr = ConstraintBuilder(pdb_file, gaps, metal)
	initial_res, final_res = constr.parse_atoms()
	constraints = constr.build_constraint(initial_res, final_res)
	return constraints




