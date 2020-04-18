import argparse
import os
import re
import time
from pele_platform.Utilities.Helpers import yaml_parser as yl

def prep_complex(complex, input_file="input.yaml", prep_output="", final_output="", debug=False):
    prep_output = run_prepwizard(complex, prep_output, debug)
    output_processed = change_repited_chain(complex, input_file, final_output)
    return output_processed
    


def run_prepwizard(complex, prep_output="", debug=False):

    input_path = os.path.abspath(complex)
    if not prep_output:
        prep_output = os.path.basename(complex.replace(".pdb", "_prep.pdb"))
    schrodinger_path = "$SCHRODINGER/utilities/prepwizard"
    
    # Run Protein Preparation Wizard - delete waters, fill missing loops and side chains
    wizard_command = "{} -fillloops -fillsidechains -delwater_hbond_cutoff 5 {} {}".format(schrodinger_path, input_path, prep_output)
    print(wizard_command)
    if not debug:
        os.system(wizard_command)
        while not os.path.exists(prep_output):
            time.sleep(10)
    return prep_output

def change_repited_chain_and_atoms(complex, residue, chain, final_output=""):
    if not final_output:
        final_output = os.path.basename(complex.replace(".pdb", "_final.pdb"))
    input_resname = residue
    input_chain = chain

    with open(complex, "r") as pdb:
        lines = pdb.readlines()
        all_chains = [line[21:22].strip() for line in lines if line.startswith("ATOM") or line.startswith("HETATM")]
        lig_chains = [line[21:22].strip() for line in lines if line[17:20] == input_resname]
        lig_atomnames = [line[12:16].strip() for line in lines if line [17:20] == input_resname]
    
    # Change ligand chain and atom names, if not unique
    lig_occ = all_chains.count(lig_chains[0])
    lig_length = len(lig_chains)
    
    unique_chain = True if lig_occ == lig_length else False
    unique_atomname = True if len(set(lig_atomnames)) == len(lig_atomnames) else False


    with open(complex, "r") as file:
        lines = file.readlines()
        dict = {} 
        atom_map = {}
        with open(final_output, "w") as final_output:
            for line in lines:
                if unique_chain == False:
                    if line[17:20] == input_resname:
                        repl_chain = "Z"
                        line = line[0:21] + line[21] + line[22:]
                else:
                    repl_chain = lig_chains[0]
                if unique_atomname == False:
                    if line[17:20] == input_resname:
                        pattern = r"[^0-9\s]"
                        letter = re.search(pattern, line[12:16])[0].strip()
                        if letter not in dict.keys():
                            dict[letter] = 1
                        else:
                            dict[letter] += 1
                        replacement = "{}{}".format(letter, dict[letter])
                        if len(replacement) != 4:
                            replacement = replacement + " "*(4 - len(replacement))
                        atom_map[line[12:16]] =  replacement
                        line = line.replace(line[12:16], replacement)
                final_output.write(line)
    return final_output.name, repl_chain, atom_map


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--complex", type=str, required=True, help="Protein-ligand complex (PDB file)")
    parser.add_argument("--input", type=str, required=True, help="YAML file to launch simulation (input.yaml)") 
    args = parser.parse_args()
    return os.path.abspath(args.complex), os.path.abspath(args.input)


if __name__ == "__main__":
    complex, input_file = parse_args()
    prep_complex(complex, input_file)
