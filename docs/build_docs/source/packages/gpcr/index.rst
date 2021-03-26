GPCR orthosteric site simulation
#########################################################################

This simulation aims to find the binding path
of a small molecule to the orthosteric site of a GPCR.

It will initially place the ligand all around the specified initial site. Then, a slight
bias will be performed against SASA and a small box containing the initial and orthosteric
sites will be set for simulation. High constraints will be placed on the CA to avoid
the collapse of the initial structure due to the lack of the membrane.

**Input** (further explained below):

    - protein-ligand PDB file
    - input.yaml specifying orthosteric and initial sites of atom

**Output** (further explained below):

    - ranked binding modes of the chosen small molecule with the desired protein

**Computational Time**: 3h

1. Complex Preparation
======================
   
Prepare the system with Protein Preparation Wizard. We would usually recommend protonating the protein (obligatory), deleting water molecules more than 5Å away from ligands
and ions as well as filling in missing loops and side chains. The system must contain the protein and the ligand, however,
the latter can be placed anywhere as its position will be automatically randomised all around the specified initial site atom.

2. Input Preparation
=====================

Prepare the input file ``input.yml``:

..  code-block:: yaml

    system: "complex.pdb" # Ligand-protein complex
    resname: "LIG" # Ligand residue name in system
    chain: "L" # Ligand chain ID
    seed: 1234
    gpcr_orth: true # Set defaults for GPCR simulation
    orthosteric_site: "A:114:CG" # Atom in the orthosteric site (chain ID:residue number:atom name)
    initial_site: "A:310:CD" # Atom in the initial site
    cpus: 50

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
