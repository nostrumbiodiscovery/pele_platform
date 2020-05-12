from Bio.PDB import PDBParser, PDBIO, Selection, NeighborSearch, Vector
import glob
import numpy as np
import os
import pele_platform.constants.constants as cs


def prepare_structure(protein_file, ligand_pdb, chain):
    
    to_remove = []

    # remove additional protein chains and all water molecules
    with open(protein_file, "r") as file:
        lines = file.readlines()

        for line in lines:
            if ((line.startswith("ATOM") or line.startswith("HETATM")) and line[21:22].strip() not in chain) or \
                    line.startswith("END") or line.startswith("TER") or line.startswith("CONECT"):
                to_remove.append(line)
            if (line.startswith("ATOM") or line.startswith("HETATM")) and line[17:20].strip() == "HOH":
                to_remove.append(line)

        protein = [line for line in lines if line not in to_remove]
    
    # read in ligand file
    ligand = []
    with open(ligand_pdb, "r") as ligand_file:
        lines = ligand_file.readlines()
        for line in lines:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                ligand.append(line)
    
    new_protein_file = os.path.basename(protein_file).replace(".pdb", "_prep.pdb")
    new_protein_file = os.path.abspath(new_protein_file)
    
    # join protein and ligand PDBs into new file
    with open(new_protein_file, "w+") as file:
        for line in protein:
            file.write(line)
        file.write("\n")
        for line in ligand:
            file.write(line)
 
    return new_protein_file


def add_water(refinement_input, chain, ligand_chain, n_waters=2):

    output = []
    refinement_input = glob.glob(refinement_input)
    n_inputs = len(refinement_input)
    water_coords = []
    resnums = []
    atomnums = []

    # get maximum residue and atom numbers
    with open(refinement_input[0], "r") as file:
        protein = file.readlines()

        for line in protein:
            if line.startswith("ATOM") or line.startswith("HETATM") or line.startswith("TER"):
                if line[23:27]:
                    resnums.append(line[23:27].strip())
                if line[7:11]:
                    atomnums.append(line[7:11].strip())
    
    resnums = [int(num) for num in resnums]
    max_resnum = max(resnums)
    water_resnums = []
    water_chain = chain  # from function args
    atomnum = max([int(num) for num in atomnums])+1
    water_atonums = []

    water = cs.water * n_waters * n_inputs

    for inp in range(n_inputs):
        for n in range(n_waters):
            O_coords = Vector([np.random.randint(0,100) for i in range(3)])
            H1_coords = O_coords + Vector(0.757, 0.586, 0.0)
            H2_coords = O_coords + Vector(-0.757, 0.586, 0.0)
            water_coords = water_coords + [list(O_coords)] + [list(H1_coords)] + [list(H2_coords)]

            max_resnum += 1 # each water must have a different residue number
            water_resnums = water_resnums+[max_resnum]*3
        max_resnum += 1
    
    water_atomnums = [atomnum+j for j in range(n_waters*3)] * n_waters
    
    # PDB lines - water
    water_output = []
    for atom, num, resnum, coord in zip(water, water_atomnums, water_resnums, water_coords):
        coord = ["{:7.4f}".format(c) for c in coord]
        coord = " ".join(coord)
        water_output.append(atom.format(num, water_chain, resnum, coord))
    
    sliced_water_output = []
    for i in range(0, len(water_output), n_waters*3):
        sliced_water_output.append(water_output[i:i+n_waters*3])
    
    # loop over minimisation inputs
    for inp, w in zip(refinement_input, sliced_water_output):
        new_protein_file = os.path.join(os.path.dirname(inp), os.path.basename(inp).replace(".pdb", "_water.pdb"))

        protein = []
        ligand = []

        # read in protein and ligand lines
        with open(inp, "r") as inp:
            lines = inp.readlines()

            for line in lines:
                if (line.startswith("ATOM") or line.startswith("HETATM")) and line[21:22].strip() == chain:
                    protein.append(line)
                if (line.startswith("ATOM") or line.startswith("HETATM")) and line[17:20].strip() == ligand_chain:
                    ligand.append(line)

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
        io.save(new_protein_file)
        output.append(new_protein_file)

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
