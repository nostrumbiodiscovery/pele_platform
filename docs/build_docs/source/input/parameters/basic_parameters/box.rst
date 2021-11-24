Perturbation site parameters
----------------------------

These are parameters to set up the region to perturb with PELE.

List of perturbation site parameters:

    1. `box_center <#box-center>`__
    2. `box_radius <#box-radius>`__
    3. `initial_site <#initial-site>`__
    4. `final_site <#final-site>`__
    5. `center_of_interface <#center-of-interface>`__

List of examples:

    - `Example 1 <#example-1>`__
    - `Example 2 <#example-2>`__


box-center
++++++++++

    - Description: The center of the perturbation box inside
      which the ligand will be perturbed.
    - Type: ``Integer``
    - Default: it depends on the package

    .. note::
       This parameter is set according to the Platform package that is chosen
       since it has a strong with the type of simulation that is pursued.
       However, if this parameter is set, it will prevail over the default
       settings of any package.

    .. note::
       Note that, currently, we only support a spherical box.

    .. seealso::
      `packages <../packages.html>`__,
      `Example 1 <#example-1>`__


box-radius
++++++++++

    - Description: The radius, in ansgtroms, of the perturbation box inside
      which the ligand will be perturbed.
    - Type: ``Integer``
    - Default: it depends on the package

    .. note::
       This parameter is set according to the Platform package that is chosen
       since it has a strong with the type of simulation that is pursued.
       However, if this parameter is set, it will prevail over the default
       settings of any package.

    .. note::
       Note that, currently, we only support a spherical box.

    .. seealso::
      `packages <../packages.html>`__,
      `Example 1 <#example-1>`__
