from Bio.PDB import PDBParser, PDBIO, Selection, NeighborSearch
import glob
import os


def prepare_structure(protein_file, ligand_pdb, chain):
    
    to_remove = []

    # remove additional protein chains and all water molecules
    with open(protein_file, "r") as file:
        print(protein_file)
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


def add_water(refinement_input, chain, ligand_chain):

    output = []
    refinement_input = glob.glob(refinement_input)
    n_waters = len(refinement_input)*2

    # hard-coded coordinates of water molecules
    water1_O = "HETATM {}  O   HOH {} {}     -25.445  -9.362   2.160  1.00  0.00           O\n"
    water1_H1 = "HETATM {}  H1  HOH {} {}     -24.805  -8.727   2.791  1.00  0.00           H\n"
    water1_H2 = "HETATM {}  H2  HOH {} {}     -26.081  -9.993   2.799  1.00  0.00           H\n"
    water2_O = "HETATM {}  O   HOH {} {}     -41.331 -14.611   0.035  1.00  0.00           O\n"
    water2_H1 = "HETATM {}  H1  HOH {} {}     -40.692 -13.976   0.666  1.00  0.00           H\n"
    water2_H2 = "HETATM {}  H2  HOH {} {}     -41.966 -15.243   0.674  1.00  0.00           H\n"

    water = [water1_O, water1_H1, water1_H2, water2_O, water2_H1, water2_H2]

    resnums = []
    atomnums = []

    # get maximum residue and atom numbers
    with open(refinement_input[0], "r") as file:
        protein = file.readlines()

        for line in protein:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                resnums.append(line[23:27].strip())
                atomnums.append(line[7:11].strip())

    water_chain = chain  # from function args
    resnums = [int(num) for num in resnums]
    water_resnums = [(max(resnums)+1)]*3+[(max(resnums)+2)]*3
    atomnum = max([int(num) for num in atomnums])
    water_atomnums = range(atomnum+1, atomnum+7)

    # PDB lines - water
    water_output = []
    for atom, num, resnum in zip(water, water_atomnums, water_resnums):
        water_output.append(atom.format(num, water_chain, resnum))

    # loop over minimisation inputs
    for inp in refinement_input:
        print("input", inp)
        new_protein_file = os.path.join(os.path.dirname(inp), os.path.basename(inp).replace(".pdb", "_water.pdb"))
        print("new protein file", new_protein_file)

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
            for line in water_output:
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
