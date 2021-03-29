Biased simulation
####################################

This simulation aims to find a binding mode containing a specific interaction between two atoms.

**Article** https://www.nature.com/articles/s41598-017-08445-5

**Input** (further explained below):

    - protein-ligand PDB file

**Output** (further explained below):

    - ranked binding modes of the chosen small molecule
      with the desired interaction

**Computational Time**: 3h

1. Complex Preparation
======================
   
Prepare the system with maestro (Protein Preparation Wizard) and output a complex.pdb. The complex.pdb must contain the protein-ligand in the desired initial conformation.
If the binding site is known, the ligand must be set as close as possible to the protein surface on that side of the protein.

2. Input Preparation
=====================

Prepare the input file ``input.yml``:

..  code-block:: yaml

    #####Normal simulation (any type)########
    system: 'docking2grid6n4b_thc.pdb' # Protein ligand PDB
    chain: 'L' # Ligand chain name
    resname: 'THC' # Ligand residue name
    seed: 12345
    # Distance to track along the simulation
    atom_dist:
    - "A:2:CA" # First atom (chain ID:residue number:atom name)
    - "B:3:CG" # Second atom
    cpus: 100
    out_in: true # Binding simulation
    initial_site: "A:577:N"
    final_site: "A:867:CB"
    ###############BIAS PART#######################
    spawning: epsilon # Apply bias
    epsilon: 0.25 # Level of bias ranging from 0 to 1
    bias_column: 7 # Column of the report starting by one to bias the results towards. (You may want to first launch a simulation with the default bias_column, then inspect the simulation report. Last, kill that simulation to launch another one with the optimized bias column value)

For more optional flags please refer to `optional flags <../../flags/index.html>`_


3. Run simulation
====================

To run the system launch the simulation with the following command:

``python -m pele_platform.main input.yml``

4. Output
=================

Raw output
+++++++++++++
Trajectory and report files for each simulation are located in ``working_folder/output``. That's where you can find
detailed information on each snapshot (PDB file, binding energy, metrics, etc.).

Selected poses
++++++++++++++++

**Clusters**

Upon completion of the simulation, all trajectories are clustered based on ligand heavy atom coordinates. Then, a cluster representative with the best binding energy (or metric of your choice) is selected.
Ranked cluster representatives can be found in:

``working_folder/results/clusters``

**Best snapshots**

In addition, top 100 structures with the best binding energy (or metric of your choice) are retrieved. This is done to ensure the clustering algorithm did not skip any valuable results. They are stored in:

``working_folder/results/top_poses``
