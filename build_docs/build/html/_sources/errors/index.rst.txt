Common errors
================

Check atoms
++++++++++++++++++++++

The following error message indicates that the atom(s) you chose to calculate the distance metric or constrain does not exist.


..  code-block:: console

    Check the atoms A:999:C given to calculate the distance metric.

To resolve the issue, you should:

- ensure the atom string has the correct format, i.e. ``"chain_id:residue_number:atom_name"``

- examine the protein structure in Schrödinger Maestro (or open the PDB file as text) and double-check, if the atom exists.

You can easily check residue numbers and atom names in the bottom panel in Maestro by hovering the mouse pointer over a specific atom. In this example, the correct atom string would be ``"A:106:OH"``.

.. image:: ../img/ppi_tutorial_1e.png
  :width: 400
  :align: center


FileNotFound
++++++++++++++++++

If PELE raises ``FileNotFoundError``, it probably means it cannot find one of the files specified in ``input.yaml`` such as system, ligand or RMSD reference. Make sure:

    - there are no spelling mistakes in file names

    - all required files are located in your working directory (or provide a relative path in your input file instead).

..  code-block:: console

    FileNotFoundError: [Errno 2] No such file or directory: '/home/anna/errors/file.pdb'


ValueError - ligand
++++++++++++++++++++

This error indicates that the software was not able to find the ligand in the PDB file. Make sure ``chain`` and ``resname`` flags
in your input file have correct values. Remember that ligands needs to have a unique chain ID!

..  code-block:: console

    ValueError: Something went wrong when extracting the ligand. Check residue&Chain on input


ValueError - water molecules
+++++++++++++++++++++++++++++++

PELE does not allow more than 4 water molecules to be perturbed, since it results in a lot of noise in simulation results. Make sure
``n_waters`` flag in the input file is set to a value between 1 and 4.

.. code-block:: console

    ValueError: Maximum 4 water molecules are allowed.

This is only relevant to ``n_waters`` flag used to automatically add water molecules to the system when running AquaPELE.