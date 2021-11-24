Constraint parameters
---------------------

These are parameters to add and modify constraints on the system to simulate with PELE.

List of constraint parameters:

    1. `constraint_level <#constraint-level>`__

List of examples:

    - `Example 1 <#example-1>`__
    - `Example 2 <#example-2>`__


constraint_level
++++++++++++++++

    - Description: Defines the level of constraining the protein backbone.
      The Platform constraints the backbone through its alpha carbon atoms (CAs).
      There are different levels of constraints:
        - ``level 0``: no constraints.
        - ``level 1``: terminal CAs constrained with a spring constant of ``5``
          kcal/mol, the rest of the CAs in the backbone with ``0.5`` kcal/mol
          at an interval of ``10``, i.e. every 10 residues.
        - ``level 2``: terminal CAs constrained at ``5`` kcal/mol, the rest
          of the CAs with ``2.5`` kcal/mol at the interval of ``8``.
        - ``level 3``: the whole backbone is constrained every ``5`` atoms
          with ``5`` kcal/mol.

    - Type: ``Integer``
    - Default: it depends on the package

    .. note::
       This parameter is set according to the Platform package that is chosen
       since it has a strong with the type of simulation that is pursued.
       We strongly suggest relying on the default
       settings for each package. However, in case of studying a system where the
       defaults are not optimal (more flexibility or rigidity required),
       this parameter can be set, and it will prevail over the default
       settings of any package.

    .. seealso::
      `packages <../packages.html>`__,
      `Example 1 <#example-1>`__,
      `Example 2 <#example-2>`__


Example 1
+++++++++

In this example we set an induced fit docking simulation with 30 computation
cores. Moreover, we reduce CAs constraints to increase the flexibility of the
system. Thus, we set the constraints level to 1.
This strategy might be useful when we want to consider the reallocation
of some regions of the backbone, such as loops, due to the perturbation
of a ligand.

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

    # Constraint parameters
    constraint_level: 1


Example 2
+++++++++

In this example we set an out --> in simulation with 50 computation
cores to study the entrance of a ligand in the cavity of a GPCR.
The best strategy to simulate membrane proteins that lack the membrane
is to add strong constraints to the protein backbone. We can set
``constraint_level`` to ``3`` for this purpose.

When using the out --> in package, we also need to set initial and final sites
in order to properly define the starting point and the region to explore
during the migration of our ligand. Check `box parameters <box.html>`__
to get further information about these two options.

..  code-block:: yaml

    # Required parameters
    system: 'system.pdb'
    chain: 'L'
    resname: 'LIG'

    # General parameters
    cpus: 50
    seed: 2021

    # Package selection
    out_in: True

    # Region selection
    initial_site: "A:40:CD"
    final_site: "A:187:ND2"

    # Constraint parameters
    constraint_level: 3
