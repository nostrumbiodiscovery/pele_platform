Pocket exploration
=====================

Introduction
--------------

Prepare your own pocket exploration simulation to identify the most promising pockets in your protein and obtain
putative binding sites of the ligand.

Check out our paper with a real-life application: `Monte Carlo simulations using PELE to identify a protein–protein inhibitor binding site and pose <https://pubs.rsc.org/en/content/articlelanding/2020/ra/d0ra01127d>`_

Inputs
++++++++

    - protein-ligand PDB file
    - YAML file with parameters

Default parameters
+++++++++++++++++++

The pocket exploration consists of two steps - global exploration to identify the pockets, followed by refinement to
optimize the binding modes.

Parameters for stage 1:

    - iterations: 100
    - pele steps: 8

Parameters for stage 2:

    - iterations: 1
    - pele steps: 250

Recommendations
+++++++++++++++++

    #. Expected computational time is around 24 h.
    #. Initial position of the small molecule is irrelevant, since it will be extracted and randomly placed all around the protein.


1. Complex Preparation
-------------------------
   
Prepare the system with Maestro (Protein Preparation Wizard): you must at least protonate the protein, but we also recommend
removing any crystallization artefacts and water molecules as well as filling the missing loops and side chains, if possible.

Ensure the ligand has:

     - unique chain ID
     - no atom names with spaces or single letters
     - any residue name except ``UNK``

2. Input Preparation
----------------------

Prepare the input file ``input.yml``:

..  code-block:: yaml

   system: 'docking2grid6n4b_thc.pdb'
   chain: 'L'
   resname: 'THC'
   seed: 12345
   #Distance to track along the simulation
   atom_dist:
   - "A:2:CA" #First atom to make the distance to
   - "B:3:CG" #Second atom to make the distance to
   site_finder: true

For more optional flags please refer to `optional flags <../../flags/index.html>`_.


3. Run simulation
-------------------

To run the system launch the simulation with the next command:

``python -m pele_platform.main input.yml``

4. Output
-------------

Best pockets ranked by ligand energy:

``working_folder/refinement_simulation/results/clusters``

Best snapshots ranked by ligand energy:

``working_folder/refinement_simulation/results/top_poses/``
