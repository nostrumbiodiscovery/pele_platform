import glob
import os
import pele_platform.Utilities.Helpers.Metals.metal_constraints as mc

def change_metal_charges(template_dir, forcefield, polarisation_factor, protein_file):

    # get names of all metals in the structure
    metals_in_structure, _ = mc.find_metals(protein_file)
    metals = list(set([atom.name for atom, _, _ in metals_in_structure]))

    # get Data and LocalData directories
    template_paths = []
    local_template_dir = template_dir
    global_template_dir = os.path.join(os.environ["PELE"], "Data")

    template_paths = find_metal_templates(metals, forcefield, global_template_dir)
    change_polarisation(template_paths, local_template_dir, polarisation_factor)


def find_metal_templates(metals, forcefield, global_template_dir):

    # find all metal template files in PELE Data
    template_paths = []
    
    for metal in metals:
        metal_template = "{}z".format(metal.lower())
        metal_template_file = glob.glob(os.path.join(global_template_dir, "Templates/{}/HeteroAtoms".format(forcefield), metal_template))
        if metal_template_file:
            template_paths.append(metal_template_file[0])
    
    return template_paths


def change_polarisation(template_paths, local_template_dir, polarisation_factor):

    # get the line with metal charge and divide the value by two
    for path in template_paths:
        with open(path, "r") as file_in:
            lines = file_in.readlines()

            for i in range(len(lines)):
                if "NBON" in lines[i]:
                    charge = lines[i+1][26:33]
                    charge_after = "{:7.6f}".format(float(charge)/polarisation_factor)
                    lines[i+1] = lines[i+1][:25] + charge_after + lines[i+1][34:]

        # create a path in DataLocal
        copied_path = os.path.join(local_template_dir, os.path.basename(path))

        # write to DataLocal
        with open(copied_path, "w+") as file_out:
            for line in lines:
                file_out.write(line)

