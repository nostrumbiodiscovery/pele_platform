Rank multiple putative poses
#####################################################################

The rescoring simulation aims to rank poses. In order to achieve this, we developed a curated side chain prediction and minimization algorithm
to relax the system. This simulation is generally used when having multiple plausible poses in order to discriminate between them.

**Input** (further explained below):

    - protein-ligand complex/es PDB file/s

**Output** (further explained below):

    - ranked binding modes

**Computational time**: 2h 

1. Complex Preparation
========================
   
Prepare the system consisting of protein and docked ligand with Schrödinger Protein Preparation Wizard. We would usually recommend protonating the protein (obligatory), deleting water molecules more than 5Å away from ligands
and ions as well as filling in missing loops and side chains.

Make sure the ligand has:

 - unique chain ID
 - unique PDB atom names with no spaces or single letters
 - any residue name except for ``UNK``

2. Input Preparation
=====================

Prepare the input file ``input.yml``:

..  code-block:: yaml

    system: 'poses*.pdb' # Protein-ligand poseS!
    chain: 'L' # Ligand chain ID
    resname: 'THC' # Ligand residue name
    seed: 12345
    # Distance between two atoms to track the simulation
    atom_dist:
    - "A:2:CA" # First atom (chain ID:residue number:atom name)
    - "B:3:CG" # Second atom
    cpus: 60
    rescoring: true # only relaxation part of the algorithm

For more optional flags please refer to `optative flags <../../documentation/index.html>`_

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

``working_folder/results/BestStructs``
