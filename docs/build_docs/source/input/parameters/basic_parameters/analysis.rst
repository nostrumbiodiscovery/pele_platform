Analysis parameters
-------------------

These are parameters to set up the analysis package of the Platform.

List of aquaPELE parameters:

    1. `analysis <#analysis>`__
    2. `only_analysis <#only-analysis>`__
    3. `bandwidth <#bandwidth>`__
    4. `max_top_clusters <#max-top-clusters>`__
    5. `top_clusters_criterion <#top-clusters-criterion>`__
    6. `cluster_representatives_criterion <#cluster-representatives-criterion>`__
    7. `max_top_poses <#max-top-poses>`__
    8. `min_population <#min-population>`__
    9. `clustering_coverage <#clustering-coverage>`__

List of examples:

    - `Example 1 <#example-1>`__
    - `Example 2 <#example-2>`__
    - `Example 3 <#example-3>`__


analysis
++++++++

    - Description: Whether to run or not the analysis at the end of the
      simulation.
    - Type: ``Boolean``
    - Default: ``True``

    .. seealso::
      `only_analysis <#only-analysis>`__,
      `Example 1 <#example-1>`__


only_analysis
+++++++++++++

    - Description: Analyze an existing PELE simulation without running a
      new one from scratch.
    - Type: ``Boolean``
    - Default: ``False``

    .. seealso::
      `analysis <#analysis>`__,
      `Example 2 <#example-2>`__


bandwidth
+++++++++

    - Description: Sets the cluster size for the clustering algorithm.
      If it is set to "auto" the software automatically will choose a
      right value
    - Type: ``Float`` or ``auto``
    - Default: ``auto``

    .. note::
       Note that ``auto`` mode will only run when using the Mean Shift
       algorithm as clustering method. Although it is the default method
       and the only one that is recommended, there are other clustering
       methods implemented in the Platform
       (check `advanced parameters <../advanced.html>`__ to get further
       information).

    .. note::
       When it is set to ``auto``, it will select the best bandwidth value
       to cover the percentage of all explored points that is set with the
       ``clustering_coverage`` parameter.


    .. seealso::
      `clustering_coverage <#clustering-coverage>`__,
      `Example 2 <#example-2>`__


max_top_clusters
++++++++++++++++

    - Description: Sets the maximum number of clusters to be selected
      as top.
    - Type: ``Integer``
    - Default: ``8``

    .. seealso::
      `Example 2 <#example-2>`__


top_clusters_criterion
++++++++++++++++++++++

    - Description: Sets the method of selecting top clusters, we can
      choose one of:
        - ``total_25_percentile`` - total energy 25th percentile
        - ``total_5_percentile`` - total energy 5th percentile
        - ``total_mean`` - total energy mean
        - ``total_min`` - total energy min
        - ``interaction_25_percentile`` - interaction energy 25th percentile
        - ``interaction_5_percentile`` - interaction energy 5th percentile
        - ``interaction_mean`` - interaction energy mean
        - ``interaction_min`` - interaction energy min
        - ``population`` - cluster population
    - Type: ``String``
    - Default: ``interaction_25_percentile``

    .. seealso::
      `Example 2 <#example-2>`__


cluster_representatives_criterion
+++++++++++++++++++++++++++++++++

    - Description: Sets method of selecting representative structures for each
      cluster, you can choose one of:
        - ``total_25_percentile`` - total energy 25th percentile
        - ``total_5_percentile`` - total energy 5th percentile
        - ``total_mean`` - total energy mean
        - ``total_min`` - total energy min
        - ``interaction_25_percentile`` - interaction energy 25th percentile
        - ``interaction_5_percentile`` - interaction energy 5th percentile
        - ``interaction_mean`` - interaction energy mean
        - ``interaction_min`` - interaction energy min
    - Type: ``String``
    - Default: ``interaction_min``

    .. seealso::
      `Example 2 <#example-2>`__


max_top_poses
+++++++++++++

    - Description: Sets the maximum number of top poses to be retrieved.
    - Type: ``Integer``
    - Default: ``100``

    .. seealso::
      `Example 2 <#example-2>`__


min_population
++++++++++++++

    - Description: Sets the minimum population that selected clusters
      must fulfil. It takes a value between 0 and 1. The default value
      of 0.01 implies that all selected clusters need to have a population
      above 1% of the total amount of sampled poses.
    - Type: ``Float``
    - Default: ``0.01``

    .. seealso::
      `Example 2 <#example-2>`__


clustering_coverage
+++++++++++++++++++

    - Description: Sets the minimum percentage of points that needs to be
      assigned to a top cluster when running mean shift clustering with
      automated ``bandwidth``. Thus, clustering bandwidth will keep
      increasing once covering the coverage percentage that is defined.
    - Type: ``Float``
    - Default: ``0.75``

    .. note::
       Note that this parameter is only used when the ``auto`` ``bandwidth``
       mode is set.

    .. seealso::
      `bandwidth <#bandwidth>`__,
      `Example 3 <#example-3>`__


Example 1
+++++++++

In this example we set an induced fit docking simulation with 30 computation
cores. Besides, we disable the analysis package so the simulation will run
but it will not be analyzed.

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

    # Analysis parameters
    analysis: False


Example 2
+++++++++

In this example we set an induced fit docking simulation with 30 computation
cores. However, instead of running the whole simulation from scratch, we
ask the analyze an existing simulation with the ``only_analysis`` option.
It is a useful feature when we want to reanalyze a previous simulation
changing some parameters, like shown below.

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

    # Analysis parameters
    only_analysis: True
    bandwidth: 8
    max_top_clusters: 12
    top_clusters_criterion: "population"
    cluster_representatives_criterion: "interaction_mean"
    max_top_poses: 20
    min_population: 0.005


Example 3
+++++++++

In this example we set an induced fit docking simulation with 30 computation
cores. For the analysis, we rely on the default bandwidth parameter, which
is ``auto``. This option finds the right clustering ``bandwidth`` for the
Mean Shift algorithm according to the ``clustering_coverage``. Thus, the
right ``bandwidth`` is selected to include inside top cluster selection, at
least, the percentage of points that is supplied with
the ``clustering_coverage`` parameter.

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

    # Analysis parameters
    only_analysis: True
    clustering_coverage: 0.60
