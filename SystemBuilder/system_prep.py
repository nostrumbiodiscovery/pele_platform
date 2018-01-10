import sys
import os
from schrodinger import structure as st


class SystemBuilder(st.StructureReader):

    def __init__(self, ligands, receptor):
        self.ligands = ligands
        self.receptor = receptor
        self.structures = st.StructureReader(self.ligands)

    def structs(self):
	"""
	   Desciption: Split each entry of a mae_file 
           retriving the residue_naem and the out file.

 	   Output:
  		self.struct_files_mae: list of output .mae files
		residues: list of residues of all lig in  mae files
	"""

	self.residues_name = []
	structures = []
        structures_mae = []
 	structures_pdb = []

        for i, structure in enumerate(self.structures):
           
	    self.residues_name.append([residue for residue in structure.residue][0].pdbres.strip())
       
	    str_name = "{}".format(structure.title)
            if str_name in structures:
                str_name = "ligand_{:d}".format(i)
	    try:
	    	structure.write(str_name+".mae")
		structure.write(str_name+".pdb")
	    except ValueError:
                str_name = "ligand_{}".format(i)
		structure.write(str_name+".mae")
                structure.write(str_name+".pdb")
	    finally:
		structures.append(str_name)
		structures_mae.append("{}.mae".format(str_name))
		structures_pdb.append("{}.pdb".format(str_name))

	self.struct_files_mae = [os.path.abspath(structure) for structure in structures_mae]
	self.struct_files_pdb = [os.path.abspath(structure) for structure in structures_pdb]

        return self.struct_files_mae, self.residues_name

    def systems(self):
	"""
	   Description: From one ligand structure build the 
	   complex with its receptor.

	   Output:
		self.struct_files_pdb: List of complexes.
	"""
	for structure in self.struct_files_pdb:
		with open(structure, "r") as pdb:
		    self.lines = pdb.readlines()
		    self.clean_lines = [line for line in self.lines if line.startswith("HETATM")]
		pdb_file = self.receptor + "TER\n" + "".join(self.clean_lines)+"END\n"	
		with open(structure, "w") as pdb_out:
		   pdb_out.write(pdb_file)
	return self.struct_files_pdb			
       
