Ligand parametrization using peleffy
=====================================

Introduction
------------

The main purpose of peleffy is to build the parameter files for PELE. Basically, PELE requires two files for each non-standard residue found in the system:

    * **IMPACT template**: a file containing the atom types and parameters of the ligand. Its job is to link each atom with the corresponding parameters using PDB atom names. Thus, PDB atom names in the input PDB file must match with the expected PDB atom names in the Impact file. This file intrinsically contains the information about the topology and connectivity of each residue.

    * **Rotamer library**: a file containing the branches that can rotate with respect to a central atomic core. Each branch consists in a set of consecutive rotatable bonds.

Besides, a third file with the Solvent parameters might be required when employing the OBC implicit solvent.

.. image:: https://martimunicoy.github.io/peleffy/_images/PELE_templates_scheme.png
    :width: 600
    :align: center

Prerequisites
---------------

Peleffy comes with the PELE Platform installation, so make sure you activate the correct python environment before
attempting the tutorial, e.g.

.. code-block:: bash

    conda activate pele_platform

You also might need Schrödinger software (any version or license type), if your ligand is not yet preprocessed.

Input files
-----------

You will need a PDB file containing the ligand, make sure it has:

    * unique PDB atom names
    * correct protonation
    * CONECT lines at the end of the file.

To ensure ligand atoms have **unique PDB atom names**:

    * Open your ligand file in Schrödinger
    * Select the ligand with a mouse click
    * Go to ``Build`` and click on ``Other edits -> Change atom properties``
    * Select ``PDB atom name`` from the drop down list and select ``Set unique PDB atom names within residues``.

If your PDB file is missing the **CONECT lines**:

    * Make sure you did not disable them in the export window of Maestro.
    * Check if your bond orders are assigned correctly, they can be adjusted with Protein Preparation Wizard tool.

Command line
-------------

The easiest way to parametrize your ligand is through the command line interface:

.. code-block:: bash

    python -m peleffy.main my_ligand.pdb

If you need to tweak any parameters, you can choose from several command line arguments:

+---------------------------------+--------------+-------------------------------------------------------+
| **Parameter**                   | **Argument** | **Functionality**                                     |
+---------------------------------+--------------+-------------------------------------------------------+
| ``-h``                          | None         | Show help message and exit                            |
|                                 |              |                                                       |
| ``--help``                      |              |                                                       |
+---------------------------------+--------------+-------------------------------------------------------+
| ``-f NAME``                     | string       | OpenForceField's forcefield name.                     |
|                                 |              |                                                       |
|                                 |              | Default is openff_unconstrained-1.2.0.offxml.         |
+---------------------------------+--------------+-------------------------------------------------------+
| ``-o PATH``                     | string       | Output path. Default is the current working directory |
|                                 |              |                                                       |
| ``--output PATH``               |              |                                                       |
+---------------------------------+--------------+-------------------------------------------------------+
| ``-r INT``                      | integer      | Rotamer library resolution in degrees. Default is 30. |
+---------------------------------+--------------+-------------------------------------------------------+
| ``--with_solvent``              | None         | Generate solvent parameters for OBC.                  |
+---------------------------------+--------------+-------------------------------------------------------+
| ``--as_datalocal``              | None         | Output will be saved following PELE's                 |
|                                 |              |                                                       |
|                                 |              | DataLocal hierarchy.                                  |
+---------------------------------+--------------+-------------------------------------------------------+
| ``-c NAME``                     | string       | The name of the method to use to compute charges,     |
|                                 |              |                                                       |
| ``--charge_method NAME``        |              | you can choose one from: gasteiger, am1bcc, OPLS.     |
+---------------------------------+--------------+-------------------------------------------------------+
| ``--charges_from_file PATH``    | string       | The path to the file with charges.                    |
+---------------------------------+--------------+-------------------------------------------------------+
| ``--chain CHAIN``               | string       | Chain ID of the molecule to parameterize.             |
+---------------------------------+--------------+-------------------------------------------------------+
| ``--include_terminal_rotamers`` | None         | Do not exclude terminal rotamers when building        |
|                                 |              |                                                       |
|                                 |              | the rotamer library.                                  |
+---------------------------------+--------------+-------------------------------------------------------+

**Example:** Parametrization of ligand contained in my_ligand.pdb file with openff-1.3.0 force field, rotamer resolution
of 10 degrees, am1bcc charge calculation method and custom output directory.

.. code-block:: bash

    python -m peleffy.main my_ligand.pdb -f openff_unconstrained-1.3.0.offxml -o "ligand_params" -r 10 -c am1bcc

API
---

If you are comfortable with python, you can make use of peleffy's API to parametrize your ligand. Follow the steps below
to get a general idea how it works and for more details, please refer to `peleffy documentation <https://martimunicoy.github.io/peleffy/index.html>`_.

1. The first step is always initializing a **Molecule object**, you can do it either from a PDB file or a string of SMILES

.. code-block:: python

    from peleffy.topology import Molecule

    molecule = Molecule("my_ligand.pdb")  # from a PDB file
    molecule = Molecule(smiles='c1ccc2cc3ccccc3cc2c1', hydrogens_are_explicit=False)  # from SMILES

2. Generating **rotamer library**

.. code-block:: python

    # Create a rotamer library
    rotamer_library = RotamerLibrary(molecule)
    rotamer_library.to_file('LIG.rot.assign')

3. **Parametrization** and creating the Impact file

a. with OpenFF

.. code-block:: python

    from peleffy.forcefield import OpenForceField
    from peleffy.topology import Topology
    from peleffy.template import Impact

    openff = OpenForceField('openff_unconstrained-1.2.0.offxml')
    parameters = openff.parameterize(molecule)
    topology = Topology(molecule, parameters)
    impact = Impact(topology)
    impact.to_file('ligz')

b. with OPLS2005

.. code-block:: python

    from peleffy.forcefield import OpenForceField
    from peleffy.topology import Topology
    from peleffy.template import Impact

    openff = OpenForceField('openff_unconstrained-1.2.0.offxml')
    parameters = openff.parameterize(molecule, charge_method='OPLS2005')
    topology = Topology(molecule, parameters)

    impact = Impact(topology)
    impact.to_file('ligz')

4. Generating **solvent parameters**

a. with OpenFF

.. code-block:: python

    from peleffy.solvent import OBC2

    solvent = OBC2(topology)  # use previously generated topology
    solvent.to_file('ligandParams.txt')

b. with OPLS2005

.. code-block:: python

    from peleffy.solvent import OPLSOBC

    solvent = OPLSOBC(topology)  # use previously generated topology
    solvent.to_file('ligandParams.txt')
