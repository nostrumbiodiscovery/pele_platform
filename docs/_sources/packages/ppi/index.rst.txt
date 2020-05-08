Prepare your own PPI simulation
####################################

This simulation aims to find a small molecule that
inhibitis the interfece between two protein domains.

**Article**: https://nostrumbiodiscovery.github.io/papers/Software/index.html#ppi-pele-monte-carlo-simulations-using-pele-to-identify-a-proteinprotein-inhibitor-binding-site-and-pose

**Input** (further explained below):

    - protein-protein pdb
    - center of the interface where to focus the simulation

**Output** (further explained below):

    - ranked binding modes of the chosen small molecule
      with the desired protein

**Computational time**: 9h

1. Complex Preparation
======================
   
Prepare the system with maestro (Protein Preparation Wizard)
and output a complex.pdb. The complex.pdb must contain protein-ligand. The ligand can be place anywhere as it will be automatically placed all around the protein interface by our automatic pipeline.

2. Input Preparation
=====================

Prepare the input file ``input.yml``:

..  code-block:: yaml

   system: '3d9u_prep.pdb' #protein-protein pdb
   protein: 'A' #Chain of the protein to keep
   ligand_pdb: "3d9u_ligand.pdb" #small molecule pdb
   chain: 'Z' #Chain of the ligand
   resname: 'LIG' #Residue of the ligand
   center_of_interface: "A:306:O" #center of the protein protein interface
   seed: 12345
   cpus: 50
   #Distance to track along the simulation
   atom_dist:
   - "A:2:CA" #First atom to make the distance to
   - "B:3:CG" #Second atom to make the distance to
   ppi: true

For more optional flags please refer to `optative falgs <../../documentation/index.html>`_


3. Run simulation
====================


To run the system launch the simulation with the next command:

``python -m pele_platform.main input.yml``

4. Output
=================

Ranked best clusters:

``working_folder/refinement_simulation/results/clusters``

Ranked best snapshots:

``working_folder/refinement_simulation/results/BestStructs/``


