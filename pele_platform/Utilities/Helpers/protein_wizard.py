import argparse
import os
import re
import time
from pele_platform.Utilities.Helpers.yaml_parser import YamlParser
import pele_platform.Checker.valid_flags as vf


def prep_complex(complex, input_file="input.yaml", prep_output="", final_output="", debug=False):

    input_path = os.path.abspath(complex)
    if not prep_output:
        prep_output = os.path.basename(complex.replace(".pdb", "_prep.pdb"))
    if not final_output:
        final_output = os.path.basename(complex.replace(".pdb", "_final.pdb"))
    schrodinger_path = "$SCHRODINGER/utilities/prepwizard"
    
    # Run Protein Preparation Wizard - delete waters, fill missing loops and side chains
    wizard_command = "{} -fillloops -fillsidechains -delwater_hbond_cutoff 5 {} {}".format(schrodinger_path, input_path, prep_output)
    if not debug:
        os.system(wizard_command)

    # Get input.yaml
    yaml = os.path.abspath(input_file)
    args = YamlParser(yaml, vf.VALID_FLAGS_PLATFORM)._parse_yaml()
    input_resname = args["resname"]

    # Read in PDB file
    while not os.path.exists(prep_output):
        time.sleep(10)

    with open(prep_output, "r") as pdb:
        lines = pdb.readlines()
        all_chains = [line[21:22].strip() for line in lines if line.startswith("ATOM") or line.startswith("HETATM")]
        lig_chains = [line[21:22].strip() for line in lines if line[17:20] == input_resname]
        lig_atomnames = [line[12:16].strip() for line in lines if line [17:20] == input_resname]
    
    # Change ligand chain and atom names, if not unique
    lig_occ = all_chains.count(lig_chains[0])
    lig_length = len(lig_chains)
    
    unique_chain = True if lig_occ == lig_length else False
    unique_atomname = True if len(set(lig_atomnames)) == len(lig_atomnames) else False


    with open(prep_output, "r") as file:
        lines = file.readlines()
        dict = {} 
        with open(final_output, "w") as final_output:
            for line in lines:
                if unique_chain == False:
                    if line[17:20] == input_resname:
                        repl_chain = "Z"
                        line = line.replace(line[21],repl_chain)
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
                        line = line.replace(line[12:16], replacement)
                final_output.write(line)

    return final_output.name


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--complex", type=str, required=True, help="Protein-ligand complex (PDB file)")
    parser.add_argument("--input", type=str, required=True, help="YAML file to launch simulation (input.yaml)") 
    args = parser.parse_args()
    return os.path.abspath(args.complex), os.path.abspath(args.input)


if __name__ == "__main__":
    complex, input_file = parse_args()
    prep_complex(complex, input_file)
