aquaPELE parameters
-------------------

These are parameters to set up the aquaPELE algorithm within a PELE simulation.

List of aquaPELE parameters:

    1. `n_waters <#n-waters>`__
    2. `waters <#waters>`__
    3. `box_water <#box-water>`__
    4. `water_radius <#water-radius>`__
    5. `water_temp <#water-temp>`__

List of examples:

    - `Example 1 <#example-1>`__
    - `Example 2 <#example-2>`__


n_waters
++++++++

    - Description: Number of new water molecules that will be included
      into the system and perturbed with aquaPELE. New water molecules
      will be contained and perturbed inside the water perturbation box
      that is defined.

    - Type: ``Integer``
    - Default: ``0``

    .. note::
       Note the difference between ``n_waters`` and ``waters``. The former
       is able to introduce new water molecules into the system before
       initializing the simulation. The latter will not add new water molecules
       but it will perturb those that are already in the input structure.

    .. seealso::
      `Example 1 <#example-1>`__


waters
++++++

    - Description: Water molecules to be perturbed with aquaPELE. They
      must belong to the input structure that is supplied to the Platform.
      Water molecules must be selected by supplying their chain and residue
      number. For example, "W:15" will select the water molecule from chain "W"
      and with residue number 15. All water molecules present in the system
      can be automatically selected with ``"all_waters"``.

    - Type: list of ``String``
    - Default: ``None``

    .. note::
       Note the difference between ``n_waters`` and ``waters``. The former
       is able to introduce new water molecules into the system before
       initializing the simulation. The latter will not add new water molecules
       but it will perturb those that are already in the input structure.

    .. seealso::
      `Example 2 <#example-2>`__


box_water
+++++++++

    - Description: Perturbation box in which water molecules will be perturbed
      with aquaPELE.

    - Type: 3-dimensional array of ``Float``
    - Default: center of mass of all water molecules to perturb

    .. note::
       Note that aquaPELE's box can only be spherical.

    .. seealso::
      `Example 1 <#example-1>`__


water_radius
++++++++++++

    - Description: Radius, in angstroms, for the perturbation box in which
      water molecules will be perturbed with aquaPELE.

    - Type: ``Float``
    - Default: ``7``

    .. note::
       Note that aquaPELE's box can only be spherical.

    .. seealso::
      `Example 1 <#example-1>`__


water_temp
++++++++++

    - Description: Temperature, in Kelvin, for the internal Metropolis criterion
      of aquaPELE. The higher it is, the easier it is to accept new water
      locations, even if they increase the energy of the system. Thus, the
      higher it is the harder it is to accept the new state of the system at the
      end of the PELE step by the outer Metropolis criterion. However, a high
      temperature promotes the sampling of water molecules.

    - Type: ``Float``
    - Default: ``5000``

    .. note::
       Note the difference between ``temperature`` and ``water_temp``. The
       former affects the global Metropolis criterion that is applied
       at the end of each PELE step and decides if we accept or reject the
       new state of the system. The latter only affects the internal
       Metropolis criterion of aquaPELE which is in charge of accepting or
       rejecting each water move.

    .. seealso::
      `temperature <pele.html#temperature>`__,
      `Example 1 <#example-1>`__


Example 1
+++++++++

In this example we set an induced fit docking simulation with 30 computation
cores. Moreover, we add 2 new water molecules inside the water box that we
define. Specifically, the spherical box we define is centered at (15, -2, 9)
and has a radius of 8 angstroms. Finally we set a temperature for the
Metropolis criterion of aquaPELE equal to 7000 Kelvin.

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

    # aquaPELE parameters
    n_waters: 2
    box_water:
      - 15
      - -2
      - 9
    water_radius: 8


Example 2
+++++++++

In this example we set an induced fit docking simulation with 30 computation
cores. Moreover, we select 3 water molecules from the input structure to
be perturbed with aquaPELE.

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

    # aquaPELE parameters
    waters:
      - "W:1"
      - "W:3"
      - "W:10"
