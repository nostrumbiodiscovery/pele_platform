Prepare your own PPI simulation
=====================================

Introduction
---------------

The PPI package aims to identify possible binding sites for a small molecule to inhibits the interface between two proteins.

Check out one of our related papers: `Monte Carlo simulations using PELE to identify a protein–protein inhibitor binding site and pose <https://pubs.rsc.org/en/content/articlelanding/2020/ra/d0ra01127d>`_.

Inputs
++++++++++

    - protein-protein PDB file
    - ligand PDB file
    - input.yaml with parameters

Default parameters
++++++++++++++++++++

The simulation consists of two stages: local exploration around the protein-protein interface and binding pose refinement.

Parameters for interface exploration:

    - iterations: 1
    - pele steps: 1000

Parameters for pose refinement:

    - iterations: 20
    - pele_steps: 12

Recommendations
++++++++++++++++++

#. We recommend tracking the distance between a ligand atom and the center of interface to aid analysis of the simulation.
#. There is no need to align the ligand, it will be automatically placed all around the protein interface by our automatic pipeline.
#. We suggest using **at least 50 CPUs**.
#. Expected computational time is around 9h.

1. Complex Preparation
-------------------------
   
Protein-protein file
++++++++++++++++++++++

The PDB file needs to be preprocessed with Maestro Protein Preparation Wizard. We usually recommend protonating the
protein (obligatory), deleting water molecules more than 5Å away from ligands and ions as well as filling in missing
loops and side chains.

Ligand file
++++++++++++++++

The ligand needs to be correctly protonated, have unique chain ID and PDB atom names. It can have any residue name except
for ``UNK``.

2. Input Preparation
------------------------

Prepare the input file ``input.yml`` as show in the template below:

.. code-block:: yaml

   system: '3d9u_prep.pdb' # protein-protein pdb
   protein: 'A' # chain of the protein to keep
   ligand_pdb: "3d9u_ligand.pdb" # ligand PDB file
   chain: 'Z' # chain ID of the ligand
   resname: 'LIG' # residue name of the ligand
   center_of_interface: "A:306:O" # center of the protein-protein interface (chain ID:residue number:atom name)
   seed: 12345
   cpus: 50
   # Distance to track along the simulation
   atom_dist:
   - "A:2:CA" # First atom to make the distance to
   - "B:3:CG" # Second atom to make the distance to
   ppi: true

For more optional flags please refer to `optional flags <../../flags/index.html>`_.


3. Run simulation
---------------------

Launch the simulation with the following command:

``python -m pele_platform.main input.yml``

4. Output
--------------

Main folders
++++++++++++++++++++++++

The working directory will contain three folders:
    - **1_interface_exploration** containing inputs and outputs of the initial induced fit simulation
    - **2_refinement_simulation** containing inputs and outputs of the rescoring simulation
    - **refinement_input** with refinement simulation inputs, i.e. best energy cluster representatives from the interface exploration.

Raw output
+++++++++++++
Trajectory and report files for each simulation are located in ``working_folder/*/output/*``. That's where you can find detailed information on each snapshot (PDB file, binding energy, metrics, etc.).

Selected poses
+++++++++++++++

Clusters
**********

Upon completion of the refinement simulation, all trajectories are clustered based on ligand heavy atom coordinates. Then, a cluster representative with the best binding energy (or metric of your choice) is selected.
Ranked cluster representatives can be found in:

``working_folder/2_refinement_simulation/results/clusters``

Best snapshots
****************

In addition, top 100 structures with the best binding energy (or metric of your choice) are retrieved. This is done to ensure the clustering algorithm did not skip any valuable results. They are stored in:

``working_folder/2_refinement_simulation/results/top_poses/``
