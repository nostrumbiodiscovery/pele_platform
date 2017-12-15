import argparse
import sys
import os
import re
from schrodinger import structure

class Lig_Prep(object):

    def __init__(self, system, lig_residue, lig_chain):
        self.system= system
        self.lig_residue = lig_residue
        self.lig_chain = lig_chain


    def extract_ligand(self):
        """
        This function returns the ligand 
        of the complex of interest
        based on residue name and chain

        :param self.system: nom del pdb
        :param ligand_resiude: ligand residue
        :param self.lig_chain: ligand chain

        :output: ligand
        """
        # receptor_with_waters_text = ""

        ligand_text = ""

        pdb_name = os.path.splitext(os.path.basename(self.system))[0]
        ligand_filename = "{}_{}".format(pdb_name, "ligand.pdb")
        # ligand_filepath = os.path.join(pdb_name, "input/{}".format(ligand_filename))

        with open(self.system, 'r') as pdb_file:
            for line in pdb_file:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    if (line[21] == self.lig_chain and line[16:21].strip(' ') == self.lig_residue):
                        ligand_text += line

        if ligand_text == "":
            raise ValueError("Something went wrong when extracting the ligand. \
            Check residue and ligand chain parameters")

        with open(ligand_filename, 'w') as ligand_file:
            ligand_file.write(ligand_text)
        return ligand_filename

    def ligand_to_mae(self, ligand_file):

        """
        Convert from pdb to mae
        """

        file, ext = ligand_file.split(".")
        st = structure.StructureReader(ligand_file).next()
        ligand_mae_file = "{}.mae".format(file)
        st.write(ligand_mae_file)
        return ligand_mae_file

def prepare_ligand(system, lig_residue, lig_chain):
    preparation = Lig_Prep(system, lig_residue, lig_chain)
    ligand_pdb = preparation.extract_ligand()
    ligand_mae = preparation.ligand_to_mae(ligand_pdb)
    return ligand_mae
