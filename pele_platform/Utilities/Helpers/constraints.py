import sys
import os
import argparse
import pele_platform.Utilities.Helpers.template_builder as tb

AMINOACIDS = ["VAL", "ASN", "GLY", "LEU", "ILE",
              "SER", "ASP", "LYS", "MET", "GLN",
              "TRP", "ARG", "ALA", "THR", "PRO",
              "PHE", "GLU", "HIS", "HIP", "TYR",
              "CYS", "HID"]

NUCLEOTIDES = ["G", "U", "A", "C"]


TER_CONSTR = 5

BACK_CONSTR = 0.5

CONSTR_ATOM = '''{{ "type": "constrainAtomToPosition", "springConstant": {0}, "equilibriumDistance": 0.0, "constrainThisAtom": "{1}:{2}:{3}" }},'''

CONSTR_DIST = '''{{ "type": "constrainAtomsDistance", "springConstant": {}, "equilibriumDistance": {}, "constrainThisAtom": "{}:{}:{}", "toThisOtherAtom": "{}:{}:{}" }},'''

CONSTR_CALPHA = '''{{ "type": "constrainAtomToPosition", "springConstant": {2}, "equilibriumDistance": 0.0, "constrainThisAtom": "{0}:{1}:_CA_" }},'''

class ConstraintBuilder(object):

    def __init__(self, pdb, gaps):
        self.pdb = pdb
        self.gaps = gaps

    def parse_atoms(self, interval=10):
        residues = {}
        initial_res = None
        with open(self.pdb, "r") as pdb:
            for line in pdb:
                resname = line[16:21].strip()
                atomtype = line[11:16].strip()
                resnum = line[22:26].strip()
                chain = line[20:22].strip()
                if line.startswith("ATOM") and ((resname in AMINOACIDS and atomtype == "CA") or resname in NUCLEOTIDES):
                    try:
                        if not initial_res:
                            residues["initial"] = [chain, line[22:26].strip()]
                            initial_res = True
                            continue
                        # Apply constraint every 10 residues
                        elif int(resnum) % interval != 1:
                            residues["terminal"] = [chain, line[22:26].strip()]
                        elif int(resnum) % interval == 1 and line.startswith("ATOM") and resname in AMINOACIDS and atomtype == "CA":
                            residues[resnum] = chain
                    except ValueError:
                        continue
        return residues

    def build_constraint(self, residues, BACK_CONSTR=BACK_CONSTR, TER_CONSTR=TER_CONSTR):

        init_constr = ['''"constraints":[''', ]

        back_constr = [CONSTR_CALPHA.format(chain, resnum, BACK_CONSTR) for resnum, chain in residues.items() if resnum.isdigit()]

        gaps_constr = self.gaps_constraints()

        terminal_constr = [CONSTR_CALPHA.format(residues["initial"][0], residues["initial"][1], TER_CONSTR), CONSTR_CALPHA.format(residues["terminal"][0], residues["terminal"][1], TER_CONSTR).strip(",")]

        final_constr = ["],"]

        constraints = init_constr + back_constr + gaps_constr + terminal_constr + final_constr

        return constraints

    def gaps_constraints(self):
        #self.gaps = {}
        gaps_constr = []
        for chain, residues in self.gaps.items():
            gaps_constr = [CONSTR_ATOM.format(TER_CONSTR, chain, terminal, "_CA_") for terminals in residues for terminal in terminals]
        return gaps_constr

def retrieve_constraints(pdb_file, gaps, metal, back_constr=BACK_CONSTR, ter_constr=TER_CONSTR, interval=10):
    constr = ConstraintBuilder(pdb_file, gaps)
    residues = constr.parse_atoms(interval=interval)
    constraints = constr.build_constraint(residues, back_constr, ter_constr)
    return constraints
