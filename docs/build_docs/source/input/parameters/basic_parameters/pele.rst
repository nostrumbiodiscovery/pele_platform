PELE general parameters
-----------------------

These are parameters that affect the general behaviour of PELE.

List of PELE general parameters:

    1. `steps <#steps>`__
    2. `minimum_steps <#minimum-steps>`__
    3. `equilibration <#equilibration>`__
    4. `equilibration_steps <#equilibration-steps>`__
    5. `temperature <#temperature>`__
    6. `anm_freq <#anm-freq>`__
    7. `sidechain_freq <#sidechain-freq>`__
    8. `min_freq <#min-freq>`__
    9. `water_freq <#water-freq>`__
    10. `conformation_freq <#conformation-freq>`__
    11. `forcefield <#forcefield>`__
    12. `solvent <#solvent>`__

List of examples:

    - `Example 1 <#example-1>`__
    - `Example 2 <#example-2>`__
    - `Example 3 <#example-3>`__
    - `Example 4 <#example-4>`__
    - `Example 5 <#example-5>`__
    - `Example 6 <#example-6>`__
    - `Example 7 <#example-7>`__

.. warning::
   Note that these parameters will not affect fragPELE.


steps
+++++

    - Description: The number of PELE steps to perform in each PELE iteration.
      Thus, this is the exact number of steps that each PELE explorer will
      execute and, as a result, all trajectories will have the same length.
    - Type: ``Integer``
    - Default: it depends on the package

    .. note::
       This parameter is set according to the Platform package that is chosen
       since it has a strong connection with the type of simulation that is pursued.
       However, if this parameter is set, it will prevail over the default
       settings of any package.

    .. note::
       Usually, the right number of PELE steps also depends on the number of
       Adaptive iterations that are set.

    .. seealso::
      `iterations <adaptive.html#iterations>`__,
      `minimum_steps <#minimum-steps>`__,
      `packages <../packages.html>`__,
      `PELE basis <../../../pele/index.html>`__,
      `Example 1 <#example-1>`__,
      `Example 2 <#example-2>`__


minimum_steps
+++++++++++++

    - Description: When this parameter is set to True, each PELE explorer
      will continue running PELE steps until all of them have reached
      the requested number of steps. As a result, those explorers that
      are faster might end up performing more steps than the minimum
      threshold.
    - Type: ``Boolean``
    - Default: ``False``

    .. note::
       This strategy allows for more efficient exploration, since the cores that already reached the
       required number of steps do not wait idly but continue the computation until all explorers have
       finished. However, it will produce an unbalanced sampling since some explorers might produce larger
       trajectories due to the fact that we force a minimum of steps to be executed but not a maximum.

    .. seealso::
      `steps <#steps>`__,
      `PELE basis <../../../pele/index.html>`__,
      `Example 2 <#example-2>`__


equilibration
+++++++++++++

    - Description: When set to true, it will equilibrate the system
      and **generate multiple starting poses**. This strategy reduces
      the bias towards the initial ligand position since the production
      run starts from the different poses that are obtained during the
      equilibration.
    - Type: ``Boolean``
    - Default: ``False``

    .. note::
       Do not confuse equilibration with pre-equilibration. The former entails
       running several equilibration steps to produce different initial
       structures. The latter only checks the amount of contacts between the
       ligand and the protein to correctly set the right clustering conditions
       for Adaptive.

    .. seealso::
      `equilibration_steps <#equilibration-steps>`__,
      `clustering_conditions <adaptive.html#clustering-conditions>`__,
      `Example 1 <#example-1>`__


equilibration_steps
+++++++++++++++++++

    - Description: The number of PELE steps to perform during the equilibration
      stage.
    - Type: ``Integer``
    - Default: ``2``

    .. note::
       This parameter will only be effective if equilibration is activated.

    .. seealso::
      `equilibration <#equilibration>`__,
      `Example 1 <#example-1>`__

temperature
+++++++++++

    - Description: The temperature in Kelvin to be used in the Metropolis
      criterion of PELE.
    - Type: ``Float``
    - Default: ``1500``

    .. seealso::
      `PELE basis <../../../pele/index.html>`__,
      `Example 3 <#example-3>`__


anm_freq
++++++++

    - Description: The frequency for the ANM algorithm of PELE. For example,
      a frequency of 1 means that it will run at every PELE step, and a
      frequency of 2 means running every 2 steps. Thus, increasing the
      frequency of the ANM algorithm will reduce the protein perturbation
      but the simulation will run faster.
    - Type: ``Integer``
    - Default: it depends on the package

    .. note::
       This parameter is set according to the Platform package that is chosen
       since it has a strong connection with the type of simulation that is pursued.
       However, if this parameter is set, it will prevail over the default
       settings of any package.

    .. note::
       Setting a frequency of 0 completely disables the ANM algorithm.

    .. seealso::
      `PELE basis <../../../pele/index.html>`__,
      `Example 3 <#example-3>`__


sidechain_freq
++++++++++++++

    - Description: The frequency for the side chain prediction algorithm
      of PELE. For example, a frequency of 1 means that it will run at every
      PELE step, and a frequency of 2 means running every 2 steps.
      Thus, increasing the frequency of the side chain prediction algorithm
      will reduce the side chain relaxation but the simulation will run
      faster.
    - Type: ``Integer``
    - Default: it depends on the package

    .. note::
       This parameter is set according to the Platform package that is chosen
       since it has a strong connection with the type of simulation that is pursued.
       However, if this parameter is set, it will prevail over the default
       settings of any package.

    .. note::
       Setting a frequency of 0 completely disables the side chain prediction
       algorithm.

    .. seealso::
      `PELE basis <../../../pele/index.html>`__,
      `Example 3 <#example-3>`__


min_freq
++++++++

    - Description: The frequency for the minimization algorithm of PELE. For example,
      a frequency of 1 means that it will run at every PELE step, and a
      frequency of 2 means running every 2 steps. Thus, increasing the
      frequency of the minimization algorithm will reduce the acceptance ratio
      of the Metropolis criterion but the simulation will run faster.
    - Type: ``Integer``
    - Default: it depends on the package

    .. note::
       This parameter is set according to the Platform package that is chosen
       since it has a strong connection with the type of simulation that is pursued.
       However, if this parameter is set, it will prevail over the default
       settings of any package.

    .. note::
       Setting a frequency of 0 completely disables the minimization algorithm.

    .. seealso::
      `PELE basis <../../../pele/index.html>`__,
      `Example 3 <#example-3>`__


water_freq
++++++++++

    - Description: The frequency for the aquaPELE algorithm of PELE. For example,
      a frequency of 1 means that it will run at every PELE step, and a
      frequency of 2 means running every 2 steps. Thus, increasing the
      frequency of aquaPELE algorithm will reduce water sampling
      but the simulation will run faster.
    - Type: ``Integer``
    - Default: it depends on the package

    .. note::
       This parameter is set according to the Platform package that is chosen
       since it has a strong connection with the type of simulation that is pursued.
       However, if this parameter is set, it will prevail over the default
       settings of any package.

    .. note::
       Note that aquaPELE is enabled only when we set some water molecules
       to be perturbed. Refer to `water parameters <water.html>`__ in order
       to get further information about how to set up aquaPELE.

    .. note::
       Setting a frequency of 0 completely disables the aquaPELE algorithm.

    .. seealso::
      `aquaPELE parameters <water.html>`__,
      `Example 4 <#example-4>`__


conformation_freq
+++++++++++++++++

    - Description: The frequency for the conformation perturbation algorithm
      of PELE. For example, a frequency of 1 means that it will run at every
      PELE step, and a frequency of 2 means running every 2 steps. Thus,
      increasing the frequency of conformation perturbation will promote the
      sampling of the different ligand conformations but the acceptance
      ratio of PELE steps might significantly decrease.
    - Type: ``Integer``
    - Default: ``4``

    .. note::
       Note that conformation perturbation is enabled only when we set
       the ``ligand_conformations`` parameter. Refer
       to `ligand parameters <ligand.html#ligand-conformations>`__ in order
       to get further information about how to set it up.

    .. note::
       Setting a frequency of 0 completely disables the conformation
       perturbation algorithm.

    .. seealso::
      `ligand_conformations <ligand.html#ligand-conformations>`__,
      `Example 5 <#example-5>`__


forcefield
++++++++++

    - Description: The force field to use during the PELE simulation. There
      are several options available:
        - ``OPLS2005``
        - ``openff-2.0.0``
        - ``openff-1.3.0``
        - ``openff-1.2.1``
        - ``openff-1.2.0``
        - ``openff-1.1.1``
        - ``openff-1.1.0``
        - ``openff-1.0.1``
        - ``openff-1.0.0``

    - Type: ``String``
    - Default: ``OPLS2005``

    .. warning::
       Selecting any OpenFF force field requires the use of peleffy to
       parametrize non standard residues. Currently, peleffy is not the
       default parametrization tool. To know how to enable it,
       check `parametrization <parametrization.html>`__ options.

    .. note::
       Using any OpenFF force field implies modeling protein residues with
       OPLS2005 and non standard residues with OpenFF.

    .. seealso::
      `use_peleffy <parametrization.html#use-peleffy>`__,
      `Example 6 <#example-6>`__


solvent
+++++++

    - Description: The implicit solvent to use during the PELE simulation.
      There are 2 options available:
        - ``VDGBNP``
        - ``OBC``

    - Type: ``String``
    - Default: ``VDGBNP`` when using ``OPLS2005`` forcefield,
      ``OBC`` when using any OpenFF force field

    .. warning::
       Note that the only implicit solvent compatible with OpenFF is ``OBC``.

    .. seealso::
      `forcefield <#forcefield>`__,
      `Example 7 <#example-7>`__


Example 1
+++++++++

In this example we set an induced fit docking simulation with 30 computation
cores. We then replace the default number of PELE steps of the induced fit
docking package. Instead of 12 steps we ask for 6. This will result in an
even faster simulation (twice as fast) at the expense of reducing the
exploration.

We are also enabling the equilibration. Thus, prior to the production run we will
run a few steps to obtain different starting positions of our ligand. The
number of PELE steps that will be devoted to the equilibration is set to 5,
replacing the default value of 2 equilibration steps.

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

    # PELE parameters
    steps: 6
    equilibration: True
    equilibration_steps: 5


Example 2
+++++++++

In this example we set an induced fit docking simulation with 30 computation
cores. We then replace the default number of PELE steps of the induced fit
docking package. Instead of 12 steps we ask for 5. Moreover, we activate
the ``minimum_steps`` mode which will transform the number of steps into a
minimum threshold. Thus, we are forcing all explorers to perform a minimum
of 5 steps but we will not block them once they finish the 5th step. Instead,
they will be able to continue executing more steps until all independent
explorers achieve the minimum threshold of 5. This strategy allows those
explorers that run faster to generate more steps, thereby increasing
the overall performance of PELE.

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

    # PELE parameters
    steps: 5
    minimum_steps: True


Example 3
+++++++++

In this example we set an induced fit docking simulation with 30 computation
cores. We then replace the default frequencies of the internal PELE algorithms.
Specifically, we are completely disabling the ANM algorithm, we are ensuring
that the side chain prediction runs at every PELE step and we are minimizing
the system every 2 steps. Finally, we are also changing the default
temperature of the Metropolis criterion, instead of 1500, we set 2000, so
the acceptance probability increases.

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

    # PELE parameters
    anm_freq: 0
    sidechain_freq: 1
    min_freq: 2
    temperature: 2000


Example 4
+++++++++

In this example we set an induced fit docking simulation with 30 computation
cores. We also activate aquaPELE by adding 2 new water molecules with
``n_waters`` parameter. Finally, we set the frequency at which aquaPELE
runs with ``water_freq`` option. Thus, it runs every 2 PELE steps.

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

    # PELE parameters
    water_freq: 2


Example 5
+++++++++

In this example we set an induced fit docking simulation with 30 computation
cores. We also provide a path that contains different conformations of
our ligand with the ``ligand_conformations`` parameter. This option will
activate the Conformation perturbation algorithm that adds an extra
perturbation step to visit all supplied ligand conformations during
the PELE simulation. However, to diminish the effects of the Conformation
perturbation algorithm, we reduce its frequency from a default of ``4``
to ``6``. This strategy modification will help to prevent the
Metropolis acceptance ratio from dropping too much.

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
    ligand_conformations: "LIG_conformations/"

    # PELE parameters
    conformation_freq: 6


Example 6
+++++++++

In this example we set an induced fit docking simulation with 30 computation
cores. We also select the latest OpenFF force field, which works with the OBC
solvent model. In order to use it, we need to activate peleffy. Check
`parametrization <parametrozation.html#use-peleffy>`__ section to get
further details.

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

    # Parametrization parameters
    use_peleffy: True

    # PELE parameters
    forcefield: "openff-2.0.0"


Example 7
+++++++++

In this example we set an induced fit docking simulation with 30 computation
cores. We also change the implicit solvent model. The default solvent when
using **OPLS2005** is **VDGBNP**. However, we can also use **OBC** if
we select it with the ``solvent`` parameter.

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

    # PELE parameters
    solvent: "OBC"
