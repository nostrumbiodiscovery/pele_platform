Prepare your own Pocket Exploration
#####################################################

To prepare your own pocket exploration to obtain putative binding sites of your small molecule and retrieve the most promising binding modes of yout ligand, please follow the steps below.

**Article**: https://nostrumbiodiscovery.github.io/papers/Software/index.html#ppi-pele-monte-carlo-simulations-using-pele-to-identify-a-proteinprotein-inhibitor-binding-site-and-pose

**Input** (further explained below):

    - Protein-ligand.pdb

**Output** (further explained below):

    - Most promising pockets ranked
    - Most promising binding modes ranked

**Computational time**: 24h

1. Complex Preparation
========================
   
Prepare the system with maestro (Protein Preparation Wizard)
and output a complex.pdb. The complex.pdb must contain protein-ligand. The ligand can be place anywhere as it will be automatically placed all around the protein by our automatic pipeline.

Make sure the ligand has:

 - Unique chain
 - No atomnames with spaces or single letter
 - Any residuename except UNK

2. Input Preparation
=====================

Prepare the input file ``input.yml``:

Pocket Exploration (24h)
+++++++++++++++++++++++++++

..  code-block:: yaml

   system: 'docking2grid6n4b_thc.pdb'
   chain: 'L'
   resname: 'THC'
   seed: 12345
   #Distance to track along the simulation
   atom_dist:
   - "A:2:CA" #First atom to make the distance to
   - "B:3:CG" #Second atom to make the distance to
   ppi: true

For more optional flags please refer to `optative falgs <../../documentation/index.html>`_


3. Run simulation
====================


To run the system launch the simulation with the next command:

``python -m pele_platform.main input.yml``

4. Output
=================

Best pockets ranked by ligand energy:

``working_folder/refinement_simulation/results/clusters``

Best snapshots ranked by ligand energy:

``working_folder/refinement_simulation/results/BestStructs/``


