PELE Converter
==============

Prerequisites
--------------
Before launching the Converter, ensure you have PELE Platform (version 1.6.2 or higher) installed and activate the Python
environment you normally use to run it.

To ensure everything is set up correctly, you can run the following command, which should print all available arguments.

.. code-block:: console

    python -m pele_platform.converter -h

Launching the tool
------------------

The most basic usage requires only **three mandatory arguments**:

    * the path to the trajectories to convert
    * format of the input files
    * desired output format.

However, you can influence the default behaviour using the parameters outlined in the `arguments section <#arguments>`_.

**Example 1.** The tool will convert all trajectories in the output folder, generated trajectories in XTC format will be saved in the same folder.

.. code-block:: console

    python -m pele_platform.converter IK1_Pele/output -if pdb -of xtc

**Example 2.** The outputs will be saved in a custom folder ``converted_trajectories``, additionally, the tool will verify the new trajectories and remove the original ones.

.. code-block:: console

    python -m pele_platform.converter IK1_Pele/output -if pdb -of xtc --delete --verify --output_path converted_trajectories

Arguments
---------
+-----------------------------+--------------+----------------------------------------------------------------+
| **Parameters**              | **Optional** | **Functionality**                                              |
+-----------------------------+--------------+----------------------------------------------------------------+
| ``input_path`` (positional) | No           | Path to the directory containing PELE trajectories to convert, |
|                             |              |                                                                |
|                             |              | this can be path the output folder, epoch folder or even       |
|                             |              |                                                                |
|                             |              | a specific trajectory file.                                    |
+-----------------------------+--------------+----------------------------------------------------------------+
| ``-if``                     | No           | The original trajectory file format, ``pdb`` or ``xtc``.       |
|                             |              |                                                                |
| ``--input_format``          |              |                                                                |
+-----------------------------+--------------+----------------------------------------------------------------+
| ``-of``                     | No           | The destination format, ``pdb`` or ``xtc``.                    |
|                             |              |                                                                |
| ``--output_format``         |              |                                                                |
+-----------------------------+--------------+----------------------------------------------------------------+
| ``-o``                      | Yes          | Path where converted trajectories will be saved.               |
|                             |              |                                                                |
| ``--output_path``           |              |                                                                |
+-----------------------------+--------------+----------------------------------------------------------------+
| ``-t``                      | Yes          | Path to the topology file in PDB format.                       |
|                             |              |                                                                |
| ``--topology``              |              |                                                                |
+-----------------------------+--------------+----------------------------------------------------------------+
| ``-d``                      | Yes          | Delete the original trajectory after conversion.               |
|                             |              | Default: False.                                                |
| ``--delete``                |              |                                                                |
+-----------------------------+--------------+----------------------------------------------------------------+
| ``-n``                      | Yes          | Number of processors to use. Default: 1.                       |
|                             |              |                                                                |
| ``--n_processors``          |              |                                                                |
+-----------------------------+--------------+----------------------------------------------------------------+
| ``--trajectory_name``       | Yes          | Name of PELE's trajectory files, e.g.                          |
|                             |              |                                                                |
|                             |              | ``trajectory.pdb`` or ``trajectory.xtc``.                      |
+-----------------------------+--------------+----------------------------------------------------------------+
| ``--verify``                | Yes          | Double-check the trajectories for consistency                  |
|                             |              |                                                                |
|                             |              | after conversion.*                                             |
+-----------------------------+--------------+----------------------------------------------------------------+
| ``--dont_verify``           | Yes          | Skip verification.*                                            |
+-----------------------------+--------------+----------------------------------------------------------------+

\* **Verification** involves checking if the original coordinates match the new ones. If they do, it will delete the original
trajectory, provided the delete option is enabled. Otherwise, a warning will be raised.
