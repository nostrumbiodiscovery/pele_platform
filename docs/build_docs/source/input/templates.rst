==========================
Custom template generation
==========================

Introduction
-------------

PELE requires two files for each non-standard residue found in the system:

    * **IMPACT template**: a file containing the atom types and parameters of the ligand. Its job is to link each atom with the corresponding parameters using PDB atom names. Thus, PDB atom names in the input PDB file must match with the expected PDB atom names in the Impact file. This file intrinsically contains the information about the topology and connectivity of each residue.

    * **Rotamer library**: a file containing the branches that can rotate with respect to a central atomic core. Each branch consists in a set of consecutive rotatable bonds.

Besides, a third file with the **Solvent parameters** might be required when employing the OBC implicit solvent.

.. image:: https://martimunicoy.github.io/peleffy/_images/PELE_templates_scheme.png
    :width: 450
    :align: center


File generation
---------------

These files are automatically created by the platform, however, the users can choose to generate them themselves
using `peleffy <../tutorials/peleffy.html>`_ via:

    * `command line interface <../tutorials/peleffy.html#command-line>`_

    * or using `Python API <../tutorials/peleffy.html#api>`_.

Once the files have been generated, you need to supply them to the PELE platform by setting specific flags
in ``input.yaml``.

# TODO: Add a link to docs about "rotamers" and "templates" flag.
