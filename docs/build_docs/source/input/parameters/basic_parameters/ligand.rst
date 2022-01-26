Ligand parameters
------------------

These are parameters that affect to how the ligand is treated during the
PELE simulation. A ligand is a small molecule that PELE perturbs
during perturbation stage in the PELE step. Check
the `PELE algorithm <../../../pele/index.html>`__ to get further information
about it.

List of ligand parameters:

    1. `ligand_conformations <#ligand-conformations>`__

List of examples:

    - `Example 1 <#example-1>`__


ligand_conformations
++++++++++++++++++++

    - Description: Defines the path to the directory containing ligand
      conformations in PDB format. When it is supplied, the conformation
      perturbation algorithm will be activated.
    - Type: ``String``
    - Default: ``None``

    .. warning::
       In order to properly apply these conformations during the PELE
       simulation, we must ensure that PDB atom names of the simulation input
       structure (``system`` parameter) match with those defined in the
       conformations supplied with ``ligand_conformations`` parameter.

    .. note::
       The path supplied with ``ligand_conformations`` parameter must
       contain several PDB files of the same ligand. They need to include
       only the ligand which must have exact PDB atom names and atom ordering.
       They only need to differ from the coordinates of their atoms thereby
       representing different ligand conformations.

    .. note::
       The number of conformations to include will depend on the chemistry
       and flexibility of the ligand. Typically, 20-30 conformations is
       a good choice in most cases.

    .. seealso::
      `Conformation perturbation <../../../pele/index.html#conformation-perturbation>`__,
      `system <../required.html>`__,
      `conformation_freq <pele.html#conformation-freq>`__,
      `Example 1 <#example-1>`__


Example 1
+++++++++

In this example we set an induced fit docking simulation with 30 computation
cores. We also provide a path that contains different conformations of
our ligand with the ``ligand_conformations`` parameter. This option will
activate the Conformation perturbation algorithm that adds an extra
perturbation step to visit all supplied ligand conformations during
the PELE simulation.

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

    # Ligand parameters
    ligand_conformations: "LIG_conformations"
