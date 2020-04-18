import urllib.request as request
import wget
import os
import time
import pele_platform.RNA.prep as rn
import pele_platform.Utilities.Helpers.protein_wizard as pw


class PDB_Fixer():

    def __init__(self, pdb, residue, chain):
        self.pdb = pdb
        self.residue = residue
        self.chain = chain
        self.pdb_proc = self.pdb

    def get_lines(self):
        with open(self.pdb_proc, "r") as f:
            return f.readlines()

    def write(self, lines):
        with open(self.pdb_proc, "w") as f:
            f.write("".join(lines))

    def fetch(self):
        try:
            wget.download(f'https://files.rcsb.org/download/{self.pdb}', self.pdb)
        except:
            raise ValueError(f"Cannot fetch {pdb} from PDB data bank. Please manually download it and use the input flag system:")

    def remove_anisou(self):
        new_lines = [line for line in self.get_lines() if not line.startswith("ANISOU")]
        self.write(new_lines)
        

    def proteinwizard(self):
        self.pdb_proc = pw.run_prepwizard(self.pdb_proc)

    def check_single_ligand_chain_name(self):
        self.pdb_proc, self.chain, self.atoms_map = pw.change_repited_chain_and_atoms(self.pdb_proc, self.residue, self.chain)

    def fix_ter(self):
        self.pdb_proc = rn.fix_ter(self.pdb_proc)

    def fix_rna(self):
        self.pdb_proc = rn.fix_rna_pdb(self.pdb_proc, self.residue)

    



def preprocess(pdb, residue, chain, rna=False):
    pdb_obj = PDB_Fixer(pdb, residue, chain)
    pdb_obj.fetch()
    pdb_obj.proteinwizard()
    pdb_obj.remove_anisou()
    pdb_obj.check_single_ligand_chain_name()
    if rna:
        pdb_obj.fix_rna()
    return pdb_obj.pdb_proc, pdb_obj.chain, pdb_obj.atoms_map
