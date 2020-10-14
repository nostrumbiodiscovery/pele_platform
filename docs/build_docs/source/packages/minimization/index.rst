Minimize your pdb
#####################################################################

The minimization package aims at minimizing the energy of the system. 
PELE will optimize all present interactions in your input pdb by slightly
perturbing protein and ligand.

**Input** (further explained below):

    - protein-ligand complex PDB file

**Output** (further explained below):

    - ranked binding modes

**Computational time**: 30min

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

    system: 'docking2grid6n4b_thc.pdb' # Protein-ligand PDB
    chain: 'L' # Ligand chain ID
    resname: 'THC' # Ligand residue name
    seed: 12345
    cpus: 20
    rescoring: true #Minimize

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
