Refine docking poses by accounting with receptor flexibility
#####################################################################

This simulation aims to enhance docking poses by having into account
receptor flexibility. To do this, a curated side chain prediction and ANM
algorithms were developed and benchmarked against standard docking techniques.

**Article**: https://www.ncbi.nlm.nih.gov/pubmed/27545443 

**Input** (further explained below):

    - Protein-ligand.pdb

**Output** (further explained below):

    - Ranked binding modes

**Computational time**: 2h 

1. Complex Preparation
========================
   
Prepare the system with maestro (Protein Preparation Wizard) and output a complex.pdb. The complex.pdb must contain the docked protein-ligand.


Make sure the ligand has:

 - Unique chain
 - No atomnames with spaces or single letter
 - Any residuename except UNK

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
    cpus: 60
    induced_fit_fast: true #2-3h less sampling but faster
    #induced_fit_exhaustive: true #6h sim but much more sampling

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
