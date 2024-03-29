=====================
Out --> in simulation
=====================

Introduction
------------------

The binding simulation aims to find the binding path of a small molecule to a given receptor. Our software will automatically
place the ligand around the specified initial site and build the simulation box based on the coordinates of the  final_site.
Upon completion of the simulation, the user will be presented with a ranking of binding modes of the small molecule.

Input
+++++++++

- protein-ligand PDB file
- YAML file with parameters

Default parameters
+++++++++++++++++++++

- iterations: 100
- pele steps: 8


Recommendations
+++++++++++++++++

#. We recommend using **at least 50 CPUs**.
#. Expected computational time is 6h.

1. Complex Preparation
--------------------------
   
Prepare the system with Maestro (Protein Preparation Wizard) and output as ``complex.pdb``.

2. Input Preparation
------------------------

Prepare the input file ``input.yml``:

..  code-block:: yaml

    system: 'docking2grid6n4b_thc.pdb' # Protein ligand PDB
    chain: 'L' # Ligand chain ID
    resname: 'THC' # Ligand residue name
    seed: 12345
    # Distance to track along the simulation
    atom_dist:
    - "A:2:CA" # First atom (chain ID:residue number:atom name)
    - "B:3:CG" # Second atom
    cpus: 60
    out_in: true
    initial_site: "A:577:N"
    final_site: "A:867:CB"

**Note:** PELE will automatically . Then simulation will start.

For more optional flags please refer to `optional flags <../../input/yaml.html>`_.


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
