Prepare your own Pocket Exploration
#####################################################

To prepare your own pocket exploration to obtain putative binding sites of your small molecule and retrieve the most
promising binding modes of your ligand, please follow the steps below.

**Article**: https://nostrumbiodiscovery.github.io/papers/Software/index.html#ppi-pele-monte-carlo-simulations-using-pele-to-identify-a-proteinprotein-inhibitor-binding-site-and-pose

**Input** (further explained below):

    - Protein-ligand.pdb

**Output** (further explained below):

    - Most promising pockets ranked
    - Most promising binding modes ranked

**Computational time**: 24h

1. Complex Preparation
========================
   
Prepare the system containing protein-ligand complex with Protein Preparation Wizard. We would usually recommend
protonating the protein (obligatory), deleting water molecules more than 5Å away from ligands and ions as well as filling in missing loops and side chains.
The ligand can be placed anywhere as its position will be automatically randomised all around the protein by our pipeline.

Make sure the ligand has:

    - unique chain ID
    - no atom names with spaces or single letter
    - any unique residue name, except ``UNK``
    - unique PDB atom names.


2. Input Preparation
=====================

Prepare the input file ``input.yml``:

..  code-block:: yaml

   system: 'docking2grid6n4b_thc.pdb'
   chain: 'L'
   resname: 'THC'
   seed: 12345
   # To track the distance between two atoms along the simulation
   atom_dist:
   - "A:2:CA" # first atom (chain ID:residue number:atom name)
   - "B:3:CG" # second atom
   allosteric: true

For more optional flags please refer to `optative flags <../../documentation/index.html>`_


3. Run simulation
====================

To run the system launch the simulation with the following command:

``python -m pele_platform.main input.yml``

4. Output
=================

Main folders
++++++++++++++++++++++++

The working directory will contain three folders:
    - **1_global_exploration** containing inputs and outputs of the global exploration
    - **2_refinement_simulation** containing inputs and outputs of the induced fit simualtion
    - **refinement_input** with refinement simulation inputs, i.e. best energy cluster representatives from the global exploration.

Raw output
+++++++++++++
Trajectory and report files for each simulation are located in ``working_folder/*/output/*``. That's where you can find detailed information on each snapshot (PDB file, binding energy, metrics, etc.).

Selected poses
+++++++++++++++
**Clusters**

Upon completion of the refinement simulation, all trajectories are clustered based on ligand heavy atom coordinates. Then, a cluster representative with the best binding energy (or metric of your choice) is selected.
Ranked cluster representatives can be found in:

``working_folder/2_refinement_simulation/results/clusters``

**Best snapshots**

In addition, top 100 structures with the best binding energy (or metric of your choice) are retrieved. This is done to ensure the clustering algorithm did not skip any valuable results. They are stored in:

``working_folder/2_refinement_simulation/results/BestStructs/``


