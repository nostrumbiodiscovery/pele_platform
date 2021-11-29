==========================================
Interaction-restricted induced fit docking
==========================================

Introduction
---------------

The induced fit simulation with interaction restrictions allow for biased exploration, where the simulation results are
limited to those that fit the specified conditions, such as specific distance or angle between user-defined atoms.

Inputs
+++++++++

    - protein-ligand PDB file
    - input.yaml specifying angle and distance restrictions

Defaults
+++++++++

    - iterations: 1
    - pele steps: 500

Recommendations
+++++++++++++++++++

    #. We recommend using **at least 50 CPUs**.
    #. Since interaction restrictions require a nested YAML file, you might encounter issues, if your YAML has incorrect indentation. Run a test simulation first (with ``test: true`` in your YAML file) to ensure everything works correctly or use one of YAML checkers available online, e.g. `YAML Lint <http://www.yamllint.com/>`_.

1. Complex Preparation
---------------------------
Prepare the system with Protein Preparation Wizard. We would usually recommend protonating the protein (obligatory), deleting water molecules more than 5Å away from ligands
and ions as well as filling in missing loops and side chains. The system must contain the protein and the ligand, however,
the latter can be placed anywhere as its position will be automatically randomised all around the specified initial site atom.

2. Input Preparation
-----------------------

Prepare the input file ``input.yml``:

Users can define two types of conditions using the atom strings (format ``chain:resnum:atomname``, e.g. A:2:CA) to select the atoms:

- **distance**: Distance between two atoms, which can be limited to a user-defined maximum, minimum or both.

- **angle**: Angle between three atoms with a user-defined maximum, minimum or both.

..  code-block:: yaml

    system: "complex.pdb" # Ligand-protein complex
    resname: "LIG" # Ligand residue name in system
    chain: "L" # Ligand chain ID

    interaction_restrictions:
    - distance:  # distance between the two atoms will not exceed 3 A
        max: 3
      atoms:
        - "A:318:OG1"   # chain A, residue number 318, atom OG1
        - "Z:201:O3"
    - angle:  # angle between those three atoms will remain between 90 and 160 degrees
        min: 90
        max: 160
      atoms:
        - "A:318:OG1"
        - "A:318:HG1"
        - "Z:201:O3"
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
