Adaptive PELE parameters
------------------------

These are parameters that affect adaptive PELE. Thus, they affect the way
how the phase space of the system is explored.

List of adaptive PELE parameters:

    1. `iterations <#iterations>`__
    2. `spawning <#spawning>`__
    3. `adaptive_restart <#adaptive-restart>`__
    4. `bias_column <#bias_column>`__
    5. `epsilon <#epsilon>`__
    6. `spawning_condition <#spawning-condition>`__
    7. `clustering_conditions <#clustering-conditions>`__
    8. `clustering_values <#clustering-values>`__

List of examples:

    - `Example 1 <#example-1>`__
    - `Example 2 <#example-2>`__
    - `Example 3 <#example-3>`__
    - `Example 4 <#example-4>`__

.. warning::
   Note that these parameters will not affect fragPELE.


iterations
++++++++++

    - Description: Adaptive iterations to run. When they are set to 1,
      Adaptive PELE will be disabled and the simulation will run with
      no exploration bias.

    - Type: ``Integer``
    - Default: it depends on the package

    .. note::
       This parameter is set according to the Platform package that is chosen
       since it has a strong connection with the type of simulation that is pursued.
       However, if this parameter is set, it will prevail over the default
       settings of any package.

    .. note::
       Adaptive iterations are also referred to as epochs.

    .. note::
       Usually, the right number of Adaptive iterations also depends on the number of
       PELE steps that are set.

    .. seealso::
      `steps <pele.html#steps>`__,
      `packages <../packages.html>`__,
      `Example 1 <#example-1>`__


spawning
++++++++

    - Description: It defines the method to spawn new structures every time
      a new Adaptive iteration starts. There are 3 options available:
        - ``independent``: Trajectories are run independently, as in
          the original PELE. It may be useful to restart simulations.
        - ``inverselyProportional``: Distributes the processors with a weight
          that is inversely proportional to the cluster population.
        - ``epsilon``: An epsilon fraction of processors are distributed
          proportionally to the value of a metric, and the rest are
          inverselyProportional distributed. This fraction must be defined
          with an additional parameter called ``epsilon``. The metric of
          interest is defined with a parameter called ``bias_column``.
          The bias towards the metric of interest depends on the
          ``spawning_condition`` parameter.

    - Type: ``String``
    - Default: ``inverselyProportional``

    .. note::
       This parameter is set according to the Platform package that is chosen
       since it has a strong connection with the type of simulation that is pursued.
       However, if this parameter is set, it will prevail over the default
       settings of any package.

    .. seealso::
      `bias_column <#bias_column>`__,
      `epsilon <#bias_column>`__,
      `spawning_condition <#spawning-condition>`__,
      `clustering_conditions <#clustering-conditions>`__,
      `packages <../packages.html>`__,
      `Example 3 <#example-3>`__


adaptive_restart
++++++++++++++++

    - Description: When this parameter is set to True, it will try to
      continue a previous simulation that might have been interrupted.
      Thus, it will start from the last Adaptive iteration that it finds
      in the working directory until the requested number of iterations
      is achieved.

    - Type: ``Boolean``
    - Default: ``True``

    .. note::
       This parameter must not be confused with ``restart``. While ``restart``
       stands for skipping any input file preparation and directly going to
       the simulation execution, it still can start from the first Adaptive
       iteration if ``adaptive_restart`` is set to False.

    .. seealso::
      `restart <general.html#restart>`__,
      `Example 4 <#example-4>`__


bias_column
+++++++++++

    - Description: Column in PELE report files that contains the metric
      of interest for Adaptive's bias. Counter starts from 1.

    - Type: ``Integer``
    - Default: it depends on the package

    .. note::
       This parameter will only be effective if ``spawning`` is set to ``epsilon``.

    .. note::
       This parameter is set according to the Platform package that is chosen
       since it has a strong connection with the type of simulation that is pursued.
       However, if this parameter is set, it will prevail over the default
       settings of any package.

    .. seealso::
      `spawning <#spawning>`__,
      `epsilon <#epsilon>`__,
      `spawning_condition <#spawning-condition>`__,
      `packages <../packages.html>`__,
      `Example 3 <#example-3>`__


epsilon
+++++++

    - Description: The fraction of the processors that will be assigned
      according to the selected metric when ``spawning`` method is set
      to ``epsilon``. It is a value between 0 and 1. The larger, the more
      bias will be applied to the metric of interest.

    - Type: ``Float``
    - Default: it depends on the package

    .. note::
       This parameter will only be effective if ``spawning`` is set to ``epsilon``.

    .. note::
       This parameter is set according to the Platform package that is chosen
       since it has a strong connection with the type of simulation that is pursued.
       However, if this parameter is set, it will prevail over the default
       settings of any package.

    .. seealso::
      `spawning <#spawning>`__,
      `bias_column <#bias_column>`__,
      `spawning_condition <#spawning-condition>`__,
      `packages <../packages.html>`__,
      `Example 3 <#example-3>`__


spawning_condition
++++++++++++++++++

    - Description: Defines how the bias towards the metric of interest
      is applied, i.e. whether it should promote clusters that minimize or
      maximize the metric of interest. There are 2 options available:
        - ``max``
        - ``min``

    - Type: ``String``
    - Default: it depends on the package

    .. note::
       This parameter will only be effective if ``spawning`` is set to ``epsilon``.

    .. note::
       This parameter is set according to the Platform package that is chosen
       since it has a strong connection with the type of simulation that is pursued.
       However, if this parameter is set, it will prevail over the default
       settings of any package.

    .. seealso::
      `spawning <#spawning>`__,
      `bias_column <#bias_column>`__,
      `epsilon <#epsilon>`__,
      `packages <../packages.html>`__,
      `Example 3 <#example-3>`__


clustering_conditions
+++++++++++++++++++++

    - Description: Defines the clustering parameters that Adaptive will employ
      to dicretize with structural clusters the conformational space of the
      ligand. The general strategy is to set up larger clusters when the ligand
      has few contacts with the protein and reduce their size when protein-ligand
      contacts increase as we want to capture this region with more detail.
      Thus, it represents an array of contacts from high to low between the
      ligand and the protein. It is related to ``clustering_values`` and the
      length of the ``clustering_conditions`` array must be equal to the length of
      ``clustering_values`` minus one.

      This parameter can be set to ``auto`` to automatically select the
      right clustering conditions. In this case, the Platform runs a
      preliminary step called pre-equilibration to capture the protein-ligand
      contacts for each particular case.

    - Type: ``List[Float]`` or ``String``
    - Default: it depends on the package

    .. note::
       This parameter is set according to the Platform package that is chosen
       since it has a strong connection with the type of simulation that is pursued.
       However, if this parameter is set, it will prevail over the default
       settings of any package.

    .. note::
       Do not confuse equilibration with pre-equilibration. The former consists
       in running several equilibration steps to produce different initial
       structures. The latter only checks the amount of contacts between the
       ligand and the protein to correctly set the right clustering conditions
       for Adaptive.

    .. seealso::
      `clustering_values <#clustering-values>`__,
      `packages <pele.html#equilibration>`__,
      `packages <../packages.html>`__,
      `Example 1 <#example-1>`__,
      `Example 2 <#example-2>`__


clustering_values
+++++++++++++++++

    - Description: Defines the clustering parameters that Adaptive will employ
      to dicretize with structural clusters the conformational space of the
      ligand. The general strategy is to set up larger clusters when the ligand
      has few contacts with the protein and reduce their size when protein-ligand
      contacts increase as we want to capture this region with more detail.
      Thus, it represents the size of each cluster, from low to high, that
      corresponds with the conditions defined in the ``clustering_conditions``
      parameter. Higher clustering values mean larger structural clusters.

    - Type: ``List[Float]``
    - Default: it depends on the package

    .. note::
       This parameter is set according to the Platform package that is chosen
       since it has a strong connection with the type of simulation that is pursued.
       However, if this parameter is set, it will prevail over the default
       settings of any package.

    .. seealso::
      `clustering_conditions <#clustering-conditions>`__,
      `packages <../packages.html>`__
      `Example 1 <#example-1>`__,
      `Example 2 <#example-2>`__


Example 1
+++++++++

In this example we set an induced fit docking simulation with 30 computation
cores. We then replace the default number of Adaptive iterations of the induced fit
docking package. Instead of 25 iterations we ask for 10. This will result in an
even faster simulation at the expense of reducing the exploration.

On the other hand, we are also specifying custom parameters for Adaptive's
clustering. We slightly reduce the sizes of clusters with the ``clustering_values``
parameter (defaults for the induced fit fast package are ``"[2.0, 5.0, 7.0]"``).
We also set ``cluster_conditions`` to ``"auto"``, so the Platform will
run a few pre-equilibration steps to determine the best cluster conditions.

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

    # Adaptive parameters
    iterations: 10
    clustering_values: "[2.0, 4.0, 6.0]"
    cluster_conditions: "auto"


Example 2
+++++++++

In this example we set an induced fit docking simulation with 30 computation
cores. We specify custom parameters for Adaptive's
clustering. We slightly reduce the sizes of clusters with the ``clustering_values``
parameter (defaults for the induced fit package are ``"[2.0, 5.0, 7.0]"``).
We also set ``cluster_conditions`` to ``"[1.5, 0.8]"``, assuming that
our ligand is able to perform more contacts than those seen in common scenarios.
For example, these conditions might work better in cases where our ligand
is highly buried in a protein cavity. Default cluster conditions for the
induced fit package are ``"[1.0, 0.6]"``.

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

    # Adaptive parameters
    clustering_values: "[2.0, 4.0, 6.0]"
    cluster_conditions: "[1.5, 0.8]"


Example 3
+++++++++

In this example we set an out --> in simulation with 50 computation
cores. When using this package, we also need to set initial and final sites
in order to properly define the starting point and the region to explore
during the migration of our ligand.
Check `perturbation site parameters <box.html>`__ to get further information
about these two options.

Then, we replace the default Adaptive spawning method of out --> in package,
which is ``inverselyProportional``, to ``epsilon``. Thus, Adaptive will
apply a certain bias towards one metric. Specifically, the portion of bias
that will be used is ``0.20``, as defined with the ``epsilon`` parameter.
Moreover, the metric of interest to track is the one in the 7th column of
PELE's reports files, which, in this case, corresponds to the distance
between the center of mass of the ligand and the chosen final site.
When setting ``spawning_condition`` to ``min``, we ask Adaptive to apply
a bit of bias towards those structures that reduce this distance, thereby
promoting the entrance of the ligand to the cavity we specified.

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
    initial_site: "A:352:CD"
    final_site: "A:283:ND2"

    # Adaptive parameters
    bias_column: 7
    spawning: "epsilon"
    epsilon: 0.20
    spawning_condition: "min"


Example 4
+++++++++

In this example we set an induced fit docking simulation with 30 computation
cores. We are disabling Adaptive restart, so in case we apply a restart, the
simulation will start from scratch, removing any Adaptive iteration that might
have been completed in a previous run.

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

    # Adaptive parameters
    adaptive_restart: False
