Perturbation site parameters
----------------------------

These are parameters to set up the region to perturb with PELE.

List of perturbation site parameters:

    1. `box_center <#box-center>`__
    2. `box_radius <#box-radius>`__
    3. `initial_site <#initial-site>`__
    4. `final_site <#final-site>`__
    5. `center_of_interface <#center-of-interface>`__

.. warning::
   Note that some parameters should only be used to define the region
   to perturb when using specific packages. They are ``initial_site``,
   ``final_site`` and ``center_of_interface``. More information is provided
   in the corresponding sections for each parameter.

List of examples:

    - `Example 1 <#example-1>`__
    - `Example 2 <#example-2>`__
    - `Example 3 <#example-3>`__
    - `Example 4 <#example-4>`__


box-center
++++++++++

    - Description: The center of the perturbation box inside
      which the ligand will be perturbed.
    - Type: 3-dimensional array of ``Float``
    - Default: it depends on the package

    .. warning::
       This parameter will prevail over the effects of other box-related
       parameters like ``initial_site``, ``final_site`` and
       ``center_of_interface`` and its combination is usually not recommended.

    .. note::
       This parameter is set according to the Platform package that is chosen
       since it has a strong connection with the type of simulation that is pursued.
       However, if this parameter is set, it will prevail over the default
       settings of any package.

    .. note::
       Note that, currently, we only support a spherical box.

    .. note::
       When custom boxes are set in a ``site_finder`` simulation, they
       change the behaviour of the package. Thus, instead of exploring
       the entire protein surface, the search will cover the custom
       box that is defined. You can check `Example 2 <#example-2>`__
       to get more details about this strategy.

    .. seealso::
      `initial_site <#initial-site>`__,
      `final_site <#final-site>`__,
      `center_of_interface <#center-of-interface>`__,
      `packages <../packages.html>`__,
      `Example 1 <#example-1>`__,
      `Example 2 <#example-2>`__


box-radius
++++++++++

    - Description: The radius, in angstroms, of the perturbation box inside
      which the ligand will be perturbed.
    - Type: ``Integer``
    - Default: it depends on the package

    .. warning::
       This parameter will prevail over the effects of other box-related
       parameters like ``initial_site``, ``final_site`` and
       ``center_of_interface`` and its combination is usually not recommended.

    .. note::
       This parameter is set according to the Platform package that is chosen
       since it has a strong connection with the type of simulation that is pursued.
       However, if this parameter is set, it will prevail over the default
       settings of any package.

    .. note::
       Note that, currently, we only support a spherical box.

    .. note::
       When custom boxes are set in a ``site_finder`` simulation, they
       change the behaviour of the package. Thus, instead of exploring
       the entire protein surface, the search will cover the custom
       box that is defined. You can check `Example 2 <#example-2>`__
       to get more details about this strategy.

    .. seealso::
      `initial_site <#initial-site>`__,
      `final_site <#final-site>`__,
      `center_of_interface <#center-of-interface>`__,
      `packages <../packages.html>`__,
      `Example 1 <#example-1>`__,
      `Example 2 <#example-2>`__


initial_site
++++++++++++

    - Description: Selection of a protein atom that is near a good starting
      point for the simulation. The Platform will place the ligand as close
      as possible to this position.
    - Type: Atom selection, ``Character``:``Integer``:``String``
    - Default: ``None``

    .. warning::
       This parameter should only be set when running the ``out_in`` package.

    .. note::
       This parameter must be used along with ``final_site`` in order
       to let the Platform properly define the perturbation box.

    .. seealso::
      `final_site <#final-site>`__,
      `packages <../packages.html>`__,
      `Example 3 <#example-3>`__,
      `Out --> In tutorial <../../../tutorials/out_in.html>`__


final_site
++++++++++

    - Description: Selection of a protein atom that is near the cavity
      that we want to visit. The Platform will build a perturbation box
      according to ``initial_site`` and ``final_site`` accordingly.
    - Type: Atom selection, ``Character``:``Integer``:``String``
    - Default: ``None``

    .. warning::
       This parameter should only be set when running the ``out_in`` package.

    .. note::
       This parameter must be used along with ``final_site`` in order
       to let the Platform properly define the perturbation box.

    .. seealso::
      `initial_site <#initial-site>`__,
      `packages <../packages.html>`__,
      `Example 3 <#example-3>`__,
      `Out --> In tutorial <../../../tutorials/out_in.html>`__


center_of_interface
+++++++++++++++++++

    - Description: Selection of a protein atom that is centered in the
      protein interface that we want the ligand to explore.
      The Platform will generate initial structures placing the ligand
      near this region and it will build a perturbation box to
      explore it properly.
    - Type: Atom selection, ``Character``:``Integer``:``String``
    - Default: ``None``

    .. warning::
       This parameter should only be set when running the ``PPI`` package.

    .. seealso::
      `packages <../packages.html>`__,
      `Example 4 <#example-4>`__,
      `PPI tutorial <../../../tutorials/ppi.html>`__


Example 1
+++++++++

In this example we set an induced fit docking simulation with 30 computation
cores. We also set a custom perturbation box. The default box used in
an induced fit docking is to center a spherical box of 6-angstrom radius
on the center of mass of the ligand. However, we can change its center
and radius with values of our choice using ``box_radius`` and ``box_center``
parameters.

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

    # Perturbation site parameters
    box_radius: 5
    box_center:
      - 21
      - -3
      - 5


Example 2
+++++++++

In this example we set a site finder simulation with 60 computation
cores. By default, this package will place our ligand all over the
protein to promote the exploration of the whole surface. However,
in case that we are interested in exploring a particular subregion of
its surface, we can define a custom box that only covers it.
Then, the Platform will spawn structures inside that box and will
guarantee that the ligand does not exit it during the simulation.

..  code-block:: yaml

    # Required parameters
    system: 'system.pdb'
    chain: 'L'
    resname: 'LIG'

    # General parameters
    cpus: 60
    seed: 2021

    # Package selection
    site_finder: True

    # Perturbation site parameters
    box_radius: 20
    box_center: 'A:241:CA'


Example 3
+++++++++

In this example we set an out --> in docking simulation with 60 computation
cores. When using this package, we also need to set initial and final sites
in order to properly define the starting point and the region to explore
during the migration of our ligand.

..  code-block:: yaml

    # Required parameters
    system: 'system.pdb'
    chain: 'L'
    resname: 'LIG'

    # General parameters
    cpus: 60
    seed: 2021

    # Package selection
    out_in: True

    # Perturbation site parameters
    initial_site: "A:43:O"
    final_site: "A:104:CD"


Example 4
+++++++++

In this example we set an out --> in docking simulation with 60 computation
cores. When using this package, we also need to set initial and final sites
in order to properly define the starting point and the region to explore
during the migration of our ligand.

..  code-block:: yaml

    # Required parameters
    system: 'system.pdb'
    chain: 'L'
    resname: 'LIG'

    # General parameters
    cpus: 60
    seed: 2021

    # Package selection
    ppi: True

    # Perturbation site parameters
    center_of_interface: "A:206:O"
