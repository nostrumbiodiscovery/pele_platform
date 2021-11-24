======================
Sidechain perturbation
======================


Introduction
---------------

Sidechain perturbation algorithm aims to optimize the interactions between a covalently bound ligand and surrounding residues.
It is composed of two stages - the initial, general docking simulation, followed by the refinement of poses with the
lowest local non-bonding energy.

Input
++++++++

    - PDB file with covalently bound ligand
    - YAML file with parameters

Default parameters
+++++++++++++++++++

Parameters for stage 1:
    - iterations: 1
    - steps: 400
    - number of structures: 1

Parameters for stage 2:
    - iterations: 1
    - steps: 100
    - number of structures: ``number of available cpus / 6`` (cannot be changed)

Recommendations
++++++++++++++++

        #. We recommend using at least 60 CPUs for the this package.

1. Complex preparation
------------------------

Prepare the system with Maestro (Protein Preparation Wizard): you must at least protonate the protein, but we also recommend
removing any crystallization artefacts and water molecules as well as filling the missing loops and side chains, if possible.

Ensure the ligand is covalently bound to a residue and its atom names remain unchanged. Unlike in other packages, the ligand
can have the same chain ID as the rest of the protein.

2. Input preparation
-----------------------

Your ``input.yaml`` should contain the following mandatory parameters:

    - **system**: path to the PDB file of the system (prepared in the previous step)
    - **resname**: residue name of the ligand
    - **chain**: chain ID of the ligand, can be the same as the protein
    - **covalent_residue**: chain and residue number of the residue to which the covalent ligand is bound, e.g. "A:13"

Additionally, you can also control:

    - **perturbation_trials**: number of trials for covalent docking. Default=10.
    - **nonbonding_radius**: the radius for calculating non-bonding energy metric. Default=20.0.
    - **refinement_angle**: during the refinement stage, the side chain perturbation will perform a very soft perturbation of a random angle between (-refinement_angle, +refinement_angle). Default=10.

..  code-block:: yaml

    system: "6M06_prepared.pdb"
    resname: "BWF"
    chain: "A"
    covalent_residue: "A:273"
    perturbation_trials: 12
    nonbonding_radius: 5
    refinement_angle: 15

For more optional flags please refer to `optional flags <../../flags/index.html>`_..


3. Run simulation
-------------------

To run the system launch the simulation with the following command:

``python -m pele_platform.main input.yml``


4. Output
------------

Best cluster representatives based on local non-bonding energy:

``working_folder/2_refinement/results/clusters``

Best snapshots:

``working_folder/2_refinement/results/top_poses/``
