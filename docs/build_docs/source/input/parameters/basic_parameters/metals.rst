Metal parameters
------------------

These are parameters that affect the treatment of metal atoms.

List of metal parameters:

    1. `polarize_metals <#polarize-metals>`__
    2. `polarization_factor <#polarization-factor>`__

List of examples:

    - `Example 1 <#example-1>`__


polarize_metals
+++++++++++++++

    - Description: Adjust charges on the metals by dividing them by
      the ``polarization_factor`` value.
    - Type: ``Boolean``
    - Default: ``False``

    .. seealso::
      `polarization_factor <#polarization-factor>`__,
      `Example 1 <#example-1>`__


polarization_factor
+++++++++++++++++++

    - Description: Factor by which metal charges should be divided when
      ``polarize_metals`` parameter is activated.
    - Type: ``Float``
    - Default: ``2``

    .. seealso::
      `polarize_metals <#polarize-metals>`__,
      `Example 1 <#example-1>`__


Example 1
+++++++++

In this example we set an induced fit docking simulation with 30 computation
cores. We also enable the polarization of metals and we ask the Platform
to divide partial charges on metal centers by ``2.2``.

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

    # Metal parameters
    polarize_metals: True
    polarization_factor: 2.2
