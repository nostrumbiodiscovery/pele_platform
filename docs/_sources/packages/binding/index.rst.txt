Prepare your own binding simulation
####################################

This simulation aims to find the binding path
of a small molecule to a given receptor.

**Input** (further explained below):

    - protein-ligand pdb

**Output** (further explained below):

    - ranked binding modes of the chosen small molecule
      with the desired protein

**Computational Time**: 6h

1. Complex Preparation
======================
   
Prepare the system with maestro (Protein Preparation Wizard)
and output a complex.pdb. The complex.pdb must contain the protein-ligand in the desired initial conformation. If the binding site is known, the ligand must be set as close as possible to the protein surface on that side of the protein.

2. Input Preparation
=====================

Prepare the input file ``input.yml``:

..  code-block:: yaml

    system: 'docking2grid6n4b_thc.pdb' #Protein ligand pdb
    chain: 'L' #Ligand chain name
    resname: 'THC' # Ligand residue name
    seed: 12345
    #Distance to track along the simulation
    atom_dist:
    - "A:2:CA" #First atom to make the distance to
    - "B:3:CG" #Second atom to make the distance to
    cpus: 100
    out_in: true

For more optional flags please refer to `optative falgs <../../documentation/index.html>`_


3. Run simulation
====================

To run the system launch the simulation with the next command:

``python -m pele_platform.main input.yml``

4. Output
=================

Best ranked clusters:

``working_folder/results/clusters``

Best ranked poses:

``working_folder/results/BestStructs/``
