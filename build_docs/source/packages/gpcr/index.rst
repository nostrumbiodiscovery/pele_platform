GPCR orthosteric site simulation
==================================

Introduction
---------------

The GPCR packe aims to find the binding path of a small molecule to the orthosteric site of a G protein-coupled receptor.

Initially, the software will place the ligand all around the user-defined initial site and create a simulation box encompassing
the initial and orthosteric sites. The simulation will be performed with a slight bias towards ligand's SASA as well as
high constraints on all alpha carbons to avoid the collapse of the structure due to the lack of the membrane. Finally,
ranked binding modes of the small molecule will be retrieved to the user.

Inputs
+++++++++

    - protein-ligand PDB file
    - input.yaml specifying orthosteric and initial sites of atom

Defaults
++++++++++++

    - iterations: 50
    - pele steps: 8
    - epsilon: 0.25
    - constraint level: 3 (backbone is constrained every 5 carbon alpha atoms with 5 kcal/mol constant)

Recommendations
+++++++++++++++++++

    #. We suggest using **at least 50 CPUs**.
    #. The simulation will take around 3h on average.


1. Complex Preparation
---------------------------
Prepare the system with Protein Preparation Wizard. We would usually recommend protonating the protein (obligatory), deleting water molecules more than 5Å away from ligands
and ions as well as filling in missing loops and side chains. The system must contain the protein and the ligand, however,
the latter can be placed anywhere as its position will be automatically randomised all around the specified initial site atom.

2. Input Preparation
-----------------------

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

For more optional flags please refer to `optional flags <../../flags/index.html>`_.


3. Run simulation
-------------------

To run the system launch the simulation with the following command:

``python -m pele_platform.main input.yml``

4. Output
--------------

Raw output
+++++++++++++
Trajectory and report files for each simulation are located in ``working_folder/output``. That's where you can find
detailed information on each snapshot (PDB file, binding energy, metrics, etc.).

Selected poses
++++++++++++++++

Clusters
**********

Upon completion of the simulation, all trajectories are clustered based on ligand heavy atom coordinates. Then, a cluster representative with the best binding energy (or metric of your choice) is selected.
Ranked cluster representatives can be found in:

``working_folder/results/clusters``

Best snapshots
*****************

In addition, top 100 structures with the best binding energy (or metric of your choice) are retrieved. This is done to ensure the clustering algorithm did not skip any valuable results. They are stored in:

``working_folder/results/top_poses``
