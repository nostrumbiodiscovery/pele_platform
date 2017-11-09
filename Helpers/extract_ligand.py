import argparse
import sys
import os


def extract_ligand(pdb_filename, general_name, ligand_chain, executing_folder):
    """
    This function creates 4 different files in .pdb format
    :param pdb_filename: nom del pdb
    :param general_name: el nom davant de l'extensio .pdb, split del nom per '.pdb'
    :param ligand_chain: cadena del lligand
    :return:
    :param executing_folder:
    """
    # receptor_with_waters_text = ""

    receptor_text = ""
    waters_text = ""
    ligand_text = ""

    ligand_filename = general_name + "_ligand.pdb"
    with open(pdb_filename, 'r') as pdb_file:
        for line in pdb_file:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                if line[21] == ligand_chain:
                    ligand_text += line
                else:
                    if line[17:20] == "HOH":
                        receptor_text += line
                    else:
                        receptor_text += line
    if ligand_text == "":
        print("Something went wrong when extracting the ligand.")
        raise 
    elif receptor_text == "":
        print("Something went wrong when extracting the receptor.")
        return False
    # else:
    #     logging.info(" - Ligand and receptor extracted correctly.")

    with open(ligand_filename, 'w') as ligand_file:
        ligand_file.write(ligand_text)
    return ligand_filename

parser = argparse.ArgumentParser()
parser.add_argument("--pdb", type=str, required=True)
parser.add_argument("--general_name", type=str, required=True)
parser.add_argument("--ligand_chain", type=str, required=True)
parser.add_argument("--executing_folder", type=str, required=True)
args = parser.parse_args()

ligand_file = extract_ligand(args.pdb, args.general_name, args.ligand_chain, args.executing_folder)

print(ligand_file)
