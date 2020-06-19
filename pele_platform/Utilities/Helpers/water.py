from Bio.PDB import PDBParser, PDBIO, Selection, NeighborSearch, Vector
import glob
import numpy as np
import os
import pele_platform.constants.constants as cs
import pele_platform.Utilities.Helpers.helpers as hp
import pele_platform.constants.pele_params as pp


def set_water_control_file(env):
    if env.waters or env.n_waters:
        env.water_arg = hp.retrieve_all_waters(
            env.system) if env.waters == "all_waters" else env.waters  # IDs of waters to perturb
        env.parameters = env.parameters.rstrip("]\n") + pp.WATER_PARAMS

        all_waters = ""
        water_string = []

        env.water_radius = env.water_radius if env.water_radius else 6
        water_inputs = [input.strip().strip('"') for input in env.adap_ex_input.split(",")]
        water_inputs = [os.path.join(env.pele_dir, elem) for elem in water_inputs]

        for inp in water_inputs:
            com = ligand_com(inp, env.residue)
            com_format = ['{:.10f}'.format(elem) for elem in com[0]]
            env.water_center = ", ".join(com_format)

            template = '{{"watersToPerturb": {{"links": {{"ids": [{index}] }}}}, "Box": {{"radius": {radius}, "fixedCenter": [{com}], "type": "sphericalBox"}}}}'

            waters = hp.retrieve_all_waters(inp)
            waters = ", ".join(['"{}"'.format(water) for water in waters])
            all_waters = all_waters + waters + ", "
            water_string.append(template.format(index=waters, radius=env.water_radius, com=env.water_center))

        env.waters = [w.strip("'") for w in water_string]
        all_waters = all_waters.rstrip(", ")
        env.water = cs.WATER.format(all_waters, env.allow_empty_selectors, env.water_temp, env.water_trials,
                                    env.water_overlap, env.water_constr, env.waters).replace("'", "")
        env.water_energy = None  # TEMPORARY FIX
    # self.water_energy = "\n".join([ cs.WATER_ENERGY.format(water.split(":")[0]) for water in water_arg ])

    else:
        env.water_energy = None
        env.water = None
        env.water_radius = None
        env.water_center = None
        env.water = ""


def water_checker(args):
    max_water = 4
    min_water = 1

    if args.n_waters:

        if args.n_waters > max_water:
            raise ValueError("Maximum {} water molecules are allowed.".format(max_water))
        elif args.n_waters < min_water:
            raise ValueError(
                "Number of water molecules (n_waters) has to be between {} and {}".format(min_water, max_water))


def add_water(refinement_input, ligand_chain, n_waters=2, test=False):
    if test:
        np.random.seed(42)

    if n_waters < 1:
        return

    else:
        output = []
        n_inputs = len(refinement_input)
        water_coords = []
        resnums = []
        atomnums = []
        chains = []
        resnames = []

        # get maximum residue and atom numbers
        with open(refinement_input[0], "r") as file:
            protein = file.readlines()

            for line in protein:
                if line.startswith("ATOM") or line.startswith("HETATM") or line.startswith("TER"):
                    try:
                        resnums.append(line[23:27].strip())
                        atomnums.append(line[7:11].strip())
                        chains.append(line[21])
                        resnames.append(line[17:20])
                    except:
                        IndexError("Line '{}' is too short".format(line))
        lig_length = resnames.count(ligand_chain)
        resnums = [int(num) for num in resnums if num]
        max_resnum = max(resnums)
        water_resnums = []
        water_chain = chains[0]  # water chain = 1st protein chain
        atomnum = max([int(num) for num in atomnums if num]) + 1 + lig_length

        water = cs.water * n_waters * n_inputs
        
        for inp in range(n_inputs):
            for n in range(n_waters):
                O_coords = Vector([np.random.randint(0, 100) for i in range(3)])
                H1_coords = O_coords + Vector(0.757, 0.586, 0.0)
                H2_coords = O_coords + Vector(-0.757, 0.586, 0.0)
                water_coords = water_coords + [list(O_coords)] + [list(H1_coords)] + [list(H2_coords)]

                max_resnum += 1  # each water must have a different residue number
                water_resnums = water_resnums + [max_resnum] * 3
            max_resnum += 1

        water_atomnums = [atomnum + j for j in range(n_waters * 3 * n_inputs)]

        # PDB lines - water
        water_output = []

        for atom, num, resnum, coord in zip(water, water_atomnums, water_resnums, water_coords):
            coord = ["{:7.4f}".format(c) for c in coord]
            coord = " ".join(coord)
            water_output.append(atom.format(num, water_chain, resnum, coord))

        sliced_water_output = []
        for i in range(0, len(water_output), n_waters * 3):
            sliced_water_output.append(water_output[i:i + n_waters * 3])

        # loop over minimisation inputs
        for inp, w in zip(refinement_input, sliced_water_output):
            new_protein_file = inp
            protein = []
            ligand = []

            # read in protein and ligand lines
            with open(inp, "r") as inp:
                lines = inp.readlines()

                for line in lines:
                    if line.startswith("ATOM") or line.startswith("HETATM"):
                        if line[17:20].strip() == ligand_chain:
                            ligand.append(line)
                        else:
                            protein.append(line)

            # add water to PDB
            with open(new_protein_file, "w+") as file:
                for line in protein:
                    file.write(line)
                file.write("\n")
                for line in w:
                    file.write(line)
                file.write("\n")
                for line in ligand:
                    file.write(line)

            # load again with Biopython
            parser = PDBParser()
            structure = parser.get_structure("complex", new_protein_file)
            water_list = []
            protein_list = Selection.unfold_entities(structure, "A")
            temp_protein_file = os.path.join(os.path.dirname(inp.name),
                                             os.path.basename(inp.name).replace(".pdb", "_temp.pdb"))

            for res in structure.get_residues():
                if res.resname == 'HOH':
                    water_list = water_list + Selection.unfold_entities(res, "A")

            # check for water contacts
            contacts5 = []
            for w in water_list:
                contacts5 = contacts5 + NeighborSearch(protein_list).search(w.coord, 5.0, "A")
            contacts5 = [c for c in contacts5 if c not in water_list]  # exclude "self" contacts

            # translate water, if needed
            while contacts5:
                contacts5 = []
                for w in water_list:
                    x, y, z = w.coord
                    w.set_coord([x - 5, y, z])
                    contacts5 = contacts5 + NeighborSearch(protein_list).search(w.coord, 5.0, "A")
                    contacts5 = [c for c in contacts5 if c not in water_list]

            # save final output
            io = PDBIO()
            io.set_structure(structure)
            io.save(temp_protein_file)
            output.append(new_protein_file)

            new_water_lines = []
            with open(temp_protein_file, "r") as temp:
                temp_lines = temp.readlines()
                for line in temp_lines:
                    if line[17:20].strip() == "HOH":
                        line = line.replace(line[7:11], str(int(line[7:11]) + lig_length))
                        if line[12:15] == "2HW":
                            line = line + "\nTER\n"
                        new_water_lines.append(line)

            new_water_lines[-2] = new_water_lines[-2].replace("\nTER\n", "")
            
            with open(new_protein_file, "w+") as file:
                for line in protein:
                    file.write(line)
                file.write("\nTER\n")
                for line in new_water_lines:
                    file.write(line)
                file.write("\n")
                for line in ligand:
                    file.write(line)
                file.write("TER")

            os.remove(temp_protein_file)

        return output


def ligand_com(refinement_input, ligand_chain):
    parser = PDBParser()
    output = []
    refinement_input = glob.glob(refinement_input)

    for inp in refinement_input:
        structure = parser.get_structure("inp", inp)
        mass = 0.0
        com = np.zeros(3)
        for res in structure.get_residues():
            if res.resname == ligand_chain:
                for atom in res.get_atoms():
                    com = com + np.array(list(atom.get_vector())) * atom.mass
                    mass += atom.mass
                com = com / mass

        output.append(com.tolist())

    return output
