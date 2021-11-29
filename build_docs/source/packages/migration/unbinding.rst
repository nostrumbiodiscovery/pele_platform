=====================
In --> out simulation
=====================

Introduction
------------------

The unbinding simulation aims to find the unbinding path of a small molecule to a given receptor.
It explores the dissociative path of a molecule and is useful when studying the dissociative binding path, or to open up the pocket
Our software will automatically set the center of the simulation box along an exit simulation.
Upon completion of the simulation, the user will be presented with a ranking of binding modes of the small molecule.

Input
+++++++++

- protein-ligand PDB file
- YAML file with parameters

Default parameters
+++++++++++++++++++++

- iterations: 100
- pele steps: 8

Modes
+++++++
- in_out
- in_out_soft

Recommendations
+++++++++++++++++

#. We recommend using **at least 30 CPUs**.
#. Expected computational time is 2h.

1. Complex Preparation
--------------------------

Prepare the system with Maestro (Protein Preparation Wizard) and output as ``complex.pdb``.

2. Input Preparation
------------------------

Prepare the input file ``input.yml``:

..  code-block:: yaml

    system: 'complex.pdb'
    chain: 'Z' # Ligand chain ID
    resname: 'STR' # Ligand residue name
    # Distance to track along the simulation
    atom_dist:
    - "A:2:HG" # First atom (chain ID:residue number:atom name
    - "A:3:CA" # Second atom
    in_out: true
    cpus: 30

For more optional flags please refer to `optional flags <../../flags/index.html>`_.

3. Run simulation
--------------------

To run the system launch the simulation with the following command:

``python -m pele_platform.main input.yml``

4. Output
----------------

Raw output
+++++++++++++
Trajectory and report files for each simulation are located in ``working_folder/output``. That's where you can find
detailed information on each snapshot (PDB file, binding energy, metrics, etc.).

Selected poses
++++++++++++++++

Clusters
***********

Upon completion of the simulation, all trajectories are clustered based on ligand heavy atom coordinates. Then, a cluster representative with the best binding energy (or metric of your choice) is selected.
Ranked cluster representatives can be found in:

``working_folder/results/clusters``

Best snapshots
******************

In addition, top 100 structures with the best binding energy (or metric of your choice) are retrieved. This is done to ensure the clustering algorithm did not skip any valuable results. They are stored in:

``working_folder/results/top_poses``
