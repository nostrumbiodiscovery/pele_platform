"""
Created on Sun Feb 23 21:23:33 2020

@author: Emiliana D'Oria
"""

# import Bio
from Bio.PDB import PDBParser, PDBIO, NeighborSearch, Selection
from Bio.PDB.vectors import rotaxis2m, Vector
import numpy as np

parser = PDBParser()

structure = parser.get_structure('protein', 'complex_ED.pdb')
ligand = parser.get_structure('ligand', 'ligand.pdb')

structure_mass = 0.
n = 0

# calculate protein COM
com = np.zeros(3)
for atom in structure.get_atoms():
    com = com + np.array(list(atom.get_vector())) * atom.mass
    structure_mass += atom.mass
com = com / structure_mass
print("Initial protein COM = {}".format(com))

# translate the protein to the origin
for atom in structure.get_atoms():
    new_pos = np.array(list(atom.get_vector())) - com
    atom.set_coord(new_pos)

# translate the ligand to the origin
original_coords = []
for atom in ligand.get_atoms():
    ligand_origin = np.array(list(atom.get_vector())) - com
    original_coords.append(ligand_origin)
    atom.set_coord(ligand_origin)
    print("Ligand at origin coords: ", ligand_origin)

# calculate protein COM after translation, should be [0, 0, 0]
com_origin = np.zeros(3)
for atom in structure.get_atoms():
    com_origin = com_origin + np.array(list(atom.get_vector())) * atom.mass
    structure_mass += atom.mass
com_origin = com_origin / structure_mass
print("COM after translation to the origin = {}".format(com_origin))

# save translated protein to file
io = PDBIO()
io.set_structure(structure)
io.save('structure_translated.pdb')

# calculating the maximum radius of the protein from the origin
coor = []
for atom in structure.get_atoms():
    coor.append(list(atom.get_vector()))

coor = np.array(coor)
coor_max = np.amax(coor, axis=0)
print("Max coordinates: ", coor_max)

d = np.sqrt(np.sum(coor_max ** 2))
print("Receptor max sphere radius is {}".format(d))

D = np.ceil(6.0 + d)  # radius of the sphere from the origin
print("Sampling {}A spherical box centered around receptor COM".format(D))
sphere_cent = com_origin

nposes = 200
j=0
print("\nGenerating {} poses...".format(nposes))
while (j < nposes):
    n += 1
    phi = np.random.uniform(0, 2 * np.pi)
    costheta = np.random.uniform(-1, 1)
    u = np.random.uniform(0, 1)
    theta = np.arccos(costheta)

    r = D * np.cbrt(u)
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)

    # move ligand to the origin
    for atom, coord in zip(ligand.get_atoms(), original_coords):
        atom.set_coord(coord)
	print("Ligand moved to origin: ", coord)
    translation = (x, y, z)
    io.set_structure(ligand)
    io.save('before_trans{}.pdb'.format(n))
    for atom in ligand.get_atoms():
        new_pos_lig_trans = np.array(list(atom.get_vector())) - translation
        atom.set_coord(new_pos_lig_trans)
        print("New ligand after translation: ", new_pos_lig_trans)
    io.set_structure(ligand)
    io.save('after_trans{}.pdb'.format(n))

    new_ligand_COM = np.zeros(3)
    ligand_mass = 0
    for atom in ligand.get_atoms():
        new_ligand_COM = new_ligand_COM + np.array(list(atom.get_vector())) * atom.mass
        ligand_mass += atom.mass
    new_ligand_COM = new_ligand_COM / ligand_mass
    print("Translated ligand COM: ", new_ligand_COM)

    # check sampling sphere
    
    dist = np.sqrt((new_ligand_COM[0] - sphere_cent[0]) ** 2 + (new_ligand_COM[1] - sphere_cent[1]) ** 2 + (
            new_ligand_COM[2] - sphere_cent[2]) ** 2)
    print("COM NEW")
    print(sphere_cent)
    print(new_ligand_COM)
    print(dist, D)

    if dist < D:
        angles = (np.random.uniform(0, 360), np.random.uniform(0, 360), np.random.uniform(0, 360))
        print("Translated ligand is inside the sampling sphere. Saved to file.")
        # random rotation
        translation = np.zeros(3)  # blank
        for atom in ligand.get_atoms():
            rotation1 = rotaxis2m(angles[0], Vector(1, 0, 0))
            atom.transform(rotation1, translation)
            rotation2 = rotaxis2m(angles[1], Vector(0, 1, 0))
            atom.transform(rotation2, translation)
            rotation3 = rotaxis2m(angles[2], Vector(0, 0, 1))
            atom.transform(rotation3, translation)
            final_coords = np.array(list(atom.get_vector()))
        print("Ligand coordinates after random rotations: ", final_coords)

	io = PDBIO()
        io.set_structure(ligand)
        io.save('newligand{}.pdb'.format(n))
        
        # check contacts at: 5A (no contacts) and 8A (needs contacts)
        contacts5 = []
        protein_list = Selection.unfold_entities(structure, "A")
        
	for atom in ligand.get_atoms():
            center = np.array(list(atom.get_vector()))
	    contacts5 = NeighborSearch(protein_list).search(center, 5.0, "A")
            contacts5_results = []
            if contacts5:
		print("Contacts found at 5A around the ligand. Aborting.")
		contacts5_results.append(False)
            else:
                print("No contacts at 5A.")
		contacts5_results.append(True)
	
	contacts8 = []
	contacts8_results = []
        
	for atom in ligand.get_atoms():
	    center = np.array(list(atom.get_vector()))
	    contacts8 = NeighborSearch(protein_list).search(center, 8.0, "A")
            if contacts8:
	        contacts8_results.append(True)
	    else:
		contacts8_results.append(False)
	if (any(contacts8_results) and all(contacts5_results)):
	    j += 1
            print("Success!")
            io.set_structure(ligand)
            io.save('correct_ligand{}.pdb'.format(n))
        else:
            print("Aborting.")
