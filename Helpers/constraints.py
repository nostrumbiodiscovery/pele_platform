import sys

AMINOACIDS = ["VAL", "ASN", "GLY", "LEU", "ILE",
              "SER", "ASP", "LYS", "MET", "GLN",
              "TRP", "ARG", "ALA", "THR", "PRO",
              "PHE", "GLU", "HIS", "HIP", "TYR",
              "CYS", "HID"]

TER_CONSTR = 5

BACK_CONSTR = 0.5

CONSTR_ATOM = '''{{ "type": "constrainAtomToPosition", "springConstant": {0}, "equilibriumDistance": 0.0, "constrainThisAtom": "{1}:{2}:{3}" }},'''

CONSTR_DIST = '''{{ "type": "constrainAtomsDistance", "springConstant": {}, "equilibriumDistance": {}, "constrainThisAtom": "{}:{}:{}", "toThisOtherAtom": "{}:{}:{}" }},'''

CONSTR_CALPHA = '''{{ "type": "constrainAtomToPosition", "springConstant": {2}, "equilibriumDistance": 0.0, "constrainThisAtom": "{0}:{1}:_CA_" }},'''

class ConstraintBuilder(object):

    def __init__(self, pdb, gaps, metals):
        self.pdb = pdb
        self.gaps = gaps
        self.metals = metals

    def parse_atoms(self):
        residues = {}
        initial_res = None
        with open(self.pdb, "r") as pdb:
            for line in pdb:
                resname = line[16:21].strip()
                atomtype = line[11:16].strip()
                resnum = line[22:26].strip()
                chain = line[20:23].strip()
                if line.startswith("ATOM") and resname in AMINOACIDS and atomtype == "CA":
                    try:
                        if not initial_res:
                            residues["initial"] = [chain, line[22:26].strip()]
                            initial_res = True
                            continue
                        # Apply constraint every 10 residues
                        elif int(resnum) % 10 != 1:
                            residues["terminal"] = [chain, line[22:26].strip()]
                        elif int(resnum) % 10 == 1 and line.startswith("ATOM") and resname in AMINOACIDS and atomtype == "CA":
                            residues[resnum] = chain
                    except ValueError:
                        continue
        return residues

    def build_constraint(self, residues):

        init_constr = ['''"constraints":[''', ]

        back_constr = [CONSTR_CALPHA.format(chain, resnum, BACK_CONSTR) for resnum, chain in residues.iteritems() if resnum.isdigit()]

        gaps_constr = self.gaps_constraints()

        metal_constr = self.metal_constraints()

        terminal_constr = [CONSTR_CALPHA.format(residues["initial"][0], residues["initial"][1], TER_CONSTR), CONSTR_CALPHA.format(residues["terminal"][0], residues["terminal"][1], TER_CONSTR).strip(",")]

        final_constr = ["],"]

        constraints = init_constr + back_constr + gaps_constr + metal_constr + terminal_constr + final_constr

        return constraints

    def gaps_constraints(self):
        #self.gaps = {}
        gaps_constr = []
        for chain, residues in self.gaps.iteritems():
            gaps_constr = [CONSTR_ATOM.format(TER_CONSTR, chain, terminal, "_CA_") for terminals in residues for terminal in terminals]
        return gaps_constr

    def metal_constraints(self):

        metal_constr = []
        for metal, ligands in self.metals.iteritems():
            metal_name, chain, metnum = metal.split(" ")
            for ligand in ligands:
                ligand_info, bond_lenght = ligand
                resname, resnum, chain, ligname = ligand_info.split(" ")
                metal_constr.append(CONSTR_DIST.format(TER_CONSTR, bond_lenght, chain, resnum, ligname, chain, metnum, metal_name))
        return metal_constr


def retrieve_constraints(pdb_file, gaps, metal):
    constr = ConstraintBuilder(pdb_file, gaps, metal)
    residues = constr.parse_atoms()
    constraints = constr.build_constraint(residues)
    return constraints
