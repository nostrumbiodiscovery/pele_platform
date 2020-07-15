import os


def remove_water(protein_file):

    new_lines = []
    new_protein_file = os.path.basename(protein_file).replace(".pdb", "_prep.pdb")
    new_protein_file = os.path.abspath(new_protein_file)

    with open(protein_file, "r") as file:
        lines = file.readlines()

        for line in lines:
            if line[17:20] != "HOH":
                new_lines.append(line)

    with open(new_protein_file, "w") as new_file:

        for line in new_lines:
            new_file.write(line)

    return new_file.name
