Simulate big conformational moves of your receptor with PCA-PELE
#####################################################################

This simulation aims to enhance the receptor movement by extracting
a pca from a movement defined for 3 or more pdbs (trajectory's snapshots).
The best binding modes of the simulated ligand will be ranked and extracted
for user inspection.


**Input** (further explained below):

    - Protein-ligand.pdb
    - 3 identical pdb (only changing coordinates) which defines a motion

**Output** (further explained below):

    - Ranked binding modes

**Computational time**: 5h 

1. Complex Preparation
========================
   
Prepare the system with maestro (Protein Preparation Wizard) and output a complex.pdb. The complex.pdb must contain the protein-ligand in the desired initial configuration.


Make sure the ligand has:

 - Unique chain
 - No atomnames with spaces or single letter
 - Any residuename except UNK

2. Input Preparation
=====================

Prepare the input file ``input.yml``:

..  code-block:: yaml

    system: 'complex.pdb' #Protein-ligand complex
    chain: 'L' #Chain name of the ligand
    resname: 'LIG' #Residue name of the ligand
    seed: 12345 
    #Calculate PCA by giving 3 pdb snapshots
    # defining a motion
    pca_traj:
    - "snap1.pdb"
    - "snap2.pdb"
    - "snap3.pdb"
    #Do not constraint structure
    remove_constraints: true
    #Big sampling of the receptor while making binding
    spawning: independent
    out_in: true #You can change the method to any other (induced_fit, rescoring...)

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
