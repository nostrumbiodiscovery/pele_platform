PELE metrics
------------

These are parameters to manage the metrics that PELE needs to track during
the simulation.

List of metric-related parameters:

    1. `atom_dist <#atom-dist>`__
    2. `rmsd_pdb <#rmsd-pdb>`__

List of examples:

    - `Example 1 <#example-1>`__
    - `Example 2 <#example-2>`__


atom_dist
+++++++++

    - Description: Asks for the distance calculation between pairs of atoms
      during the PELE simulation. More than one pair of atoms can be
      appended to the list if more than one atomic distance must be retrieved.
    - Type: list of ``Character``:``Integer``:``String````Integer``
    - Default: ``None``

    .. seealso::
      `rmsd_pdb <#rmsd-pdb>`__,
      `Example 1 <#example-1>`__,
      `Example 2 <#example-2>`__


rmsd_pdb
++++++++

    - Description: Defines the path to a PDB structure that will be used
      to calculate the RMSD of the ligand during the PELE simulation.
    - Type: ``String``
    - Default: ``None``

    .. note::
       Note that the PDB structure that is supplied must contain the same
       ligand as the input structure of PELE. The ligand must share the same
       residue name and number, and also the same atom names as the input
       structure.

    .. seealso::
      `atom_dist <#atom-dist>`__,
      `Example 1 <#example-1>`__


Example 1
+++++++++

In this example we set an induced fit docking simulation with 30 computation
cores. We also ask for the calculation of a distance between two atoms and
the RMSD of the ligand with respect to a reference structure.

..  code-block:: yaml

    # Required parameters
    system: 'system.pdb'
    chain: 'L'
    resname: 'LIG'

    # General parameters
    cpus: 30
    seed: 2021

    # Package selection
    induced_fit_fast: True

    # PELE metrics parameters
    atom_dist:
        - "A:2:CA"
        - "B:3:CG"
    rmsd_pdb: "reference.pdb"


Example 2
+++++++++

In this example we set an induced fit docking simulation with 30 computation
cores. We also ask for the calculation of two distances between atoms
A:2:CA and B:3:CG, and also between A:5:N and B:3:CG.

..  code-block:: yaml

    # Required parameters
    system: 'system.pdb'
    chain: 'L'
    resname: 'LIG'

    # General parameters
    cpus: 30
    seed: 2021

    # Package selection
    induced_fit_fast: True

    # PELE metrics parameters
    atom_dist:
        - "A:2:CA"
        - "B:3:CG"
        - "A:5:N"
        - "B:3:CG"
