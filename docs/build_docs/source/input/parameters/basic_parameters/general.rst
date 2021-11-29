General parameters
------------------

These are general parameters that affect the global behaviour of the Platform.

List of general parameters:

    1. `cpus <#cpus>`__
    2. `test <#test>`__
    3. `debug <#debug>`__
    4. `restart <#restart>`__
    5. `seed <#seed>`__

List of examples:

    - `Example 1 <#example-1>`__
    - `Example 2 <#example-2>`__
    - `Example 3 <#example-3>`__


cpus
++++

    - Description: The number of computation cores that the Platform will use.
    - Type: ``Integer``
    - Default: ``60``

    .. note::
      When running on test mode the number of cpus is ignored and the simulation
      only uses 5 of them.

    .. seealso::
      `test <general.html#test>`__,
      `Example 1 <#example-1>`__,
      `Example 2 <#example-2>`__,
      `Example 3 <#example-3>`__


test
++++

    - Description: Run a quick test to check the simulation works (~2 min).
    - Type: ``Boolean``
    - Default: ``False``

    .. warning::
       Never use the control files generated in test mode as input for a
       production simulation because other parameters such as temperature,
       ANM and minimization are tweaked to make simulation faster.

    .. seealso::
      `cpus <general.html#cpus>`__,
      `Example 1 <#example-1>`__,
      `Example 2 <#example-2>`__


debug
+++++

    - Description: Use this flag to only create the inputs of the simulation.
      No actual simulation will run.
    - Type: ``Boolean``
    - Default: ``False``

    .. note::
      The combination of this parameter along with ``debug`` allows the user to manually modify any input file that the
      Platform generates. For example, running first in debug mode will generate input files like ``pele.conf``,
      ``adaptive.conf`` or ligand templates. The user can then modify them at their will and restart the job with
      ``restart``.

    .. seealso::
      `restart <general.html#restart>`__,
      `Example 1 <#example-1>`__


restart
+++++++

    - Description: Use restart parameter set to True to start a simulation
      from scratch (with existing input PDBs and configuration files).
    - Type: ``Boolean``
    - Default: ``False``

    .. note::
      The combination of this parameter along with ``debug`` allows the user to manually modify any input file that the
      Platform generates. For example, running first on debug mode will generate input files like ``pele.conf``,
      ``adaptive.conf`` or ligand templates. The user can then modify them at their will and restart the job with
      ``restart``.

    .. note::
       This parameter must not be confused with ``adaptive_restart``.
       While ``restart`` stands for skipping any input file preparation
       and directly going to the simulation execution, it still can start
       from the first Adaptive iteration if ``adaptive_restart`` is set to
       False.

    .. seealso::
      `debug <general.html#debug>`__,
      `adaptive_restart <adaptive.html#adaptive_restart>`__,
      `Example 1 <#example-1>`__,
      `Example 3 <#example-3>`__


seed
++++

    - Description: Seed for Platform's pseudo-random numbers generator for reproducibility.
      When no ``seed`` is set, it will be initialized to a random number. This
      random number can be consulted afterwards by checking ``adaptive.conf`` file.

    - Type: ``Integer``
    - Default: ``None``

    .. note::
      It is always a good practice to establish a fixed seed in order to guarantee
      reproducibility.

    .. seealso::
      `Example 1 <#example-1>`__


working_folder
++++++++++++++

    - Description: Directory where the simulation will run.

    - Type: ``String``
    - Default: ``LIG_Pele``, where LIG is the residue name of our ligand

    .. note::
      When ``working_folder`` is not set, the default behaviour is to never
      replace an existing folder. So, in case that ``LIG_Pele`` directory
      already exists, the ultimate ``working_folder`` will be set to
      ``LIG_Pele_1``, ``LIG_Pele_2``, and so on.

    .. seealso::
      `Example 1 <#example-1>`__,
      `Example 3 <#example-3>`__


Example 1
+++++++++

In this example we set an induced fit docking simulation with 10 computation
cores and run it in debug mode.
Moreover, test and restart modes are disabled. Finally, we also establish a
specific seed for the pseudo-random numbers generator and a custom working folder.

..  code-block:: yaml

    # General parameters
    cpus: 10
    test: False
    debug: True
    restart: False
    seed: 2021
    working_folder: "my_custom_choice"

    # Required parameters
    system: 'system.pdb'
    chain: 'L'
    resname: 'LIG'

    # Package selection
    induced_fit_fast: True


Example 2
+++++++++

In this example we set an induced fit docking simulation with 10 computation
cores and run it in test mode.
When using this mode, the number of computation cores that will be used is always
going to be 5, regardless of the number of cores requested with the ``cpus`` parameter.

..  code-block:: yaml

    # General parameters
    cpus: 10
    test: True

    # Required parameters
    system: 'system.pdb'
    chain: 'L'
    resname: 'LIG'

    # Package selection
    induced_fit_fast: True


Example 3
+++++++++

In this example we ask the induced fit docking simulation to be restarted.
Consequently, the Platform expects to find a directory previously created
with valid input files. To generate them, we need to execute the Platform
in debug mode, as shown in `Example 1 <#example-1>`__. So, in this case
the working_folder that we set it must already exist.

..  code-block:: yaml

    # General parameters
    cpus: 10
    restart: True
    working_folder: "my_custom_choice"

    # Required parameters
    system: 'system.pdb'
    chain: 'L'
    resname: 'LIG'

    # Package selection
    induced_fit_fast: True
