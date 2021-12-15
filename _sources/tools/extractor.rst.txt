PELE Pathway Extractor
======================

Prerequisites
--------------
Before launching the Pathway Extractor, ensure you have PELE Platform (version 1.6.2 or higher) installed and activate the Python
environment you normally use to run it.

To ensure everything is set up correctly, you can run the following command, which should print all available arguments.

.. code-block:: console

    python -m pele_platform.pathway_extractor -h

Launching the tool
------------------

The PELE Pathway Extractor requires **three mandatory arguments**:

    * path to the epoch (iteration) folder
    * trajectory number
    * snapshot (model) number.

To further customize the functionality, please refer to the `full list of arguments <#arguments>`_.

**Example 1.** Extracting structure from epoch 0, trajectory 1, snapshot 0.

.. code-block:: console

    python -m pele_platform.pathway_extractor IK1_Pele/output/0 1 0


**Example 2.** Performing the same extraction as in the Example 1, but this time providing a custom output directory,
PDB file name and path to the topology file (mandatory for XTC trajectories).

.. code-block:: console

    python -m pele_platform.pathway_extractor IK1_Pele/output/0 1 0 --top IK1_Pele/output/topologies/topology_0.pdb -o my_dir --name pathway1.pdb

Arguments
---------

Full list of command line arguments.

+----------------+-------------------------------------------------------------------------------------+
| **Parameters** | **Functionality**                                                                   |
+----------------+-------------------------------------------------------------------------------------+
| ``epoch``      | Path to the epoch (iteration) folder.                                               |
+----------------+-------------------------------------------------------------------------------------+
| ``trajectory`` | Trajectory number of the snapshot to extract.                                       |
+----------------+-------------------------------------------------------------------------------------+
| ``snapshot``   | Snapshot to select (in accepted steps).                                             |
+----------------+-------------------------------------------------------------------------------------+
| ``-o``         | Output directory to save the resulting pathway. It will default to current working  |
|                |                                                                                     |
|                | directory.                                                                          |
+----------------+-------------------------------------------------------------------------------------+
| ``--name``     | Name of the PDB file that will be generated.Default: "pathway.pdb".                 |
+----------------+-------------------------------------------------------------------------------------+
| ``--top``      | Path to the PDB topology file. **Mandatory when extracting from XTC directories.**  |
+----------------+-------------------------------------------------------------------------------------+
