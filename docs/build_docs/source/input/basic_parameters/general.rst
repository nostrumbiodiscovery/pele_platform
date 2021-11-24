General parameters
------------------

These are general parameters that affect to the global behaviour of the Platform.


test
++++

    - Description: Run a quick test to check the simulation works (~2 min).
    - Type: ``Boolean``
    - Default: ``False``

    .. warning::
       Never use the control files from the test as input for a production simulation
       as temperature, ANM and minimization are twicked to make simulation faster.


debug
+++++

    - Description: Use this flag to only create the inputs of the simulation. No simulation is run.
    - Type: ``Boolean``
    - Default: ``False``

    .. note::
      The combination of this parameter along with ``debug`` allows the user to manually modify any input file that the
      Platform generates. For example, running first on debug mode will generate input files like ``pele.conf``,
      ``adaptive.conf`` or ligand templates. The user can then modify them at their will and restart the job with
      ``restart``.

    .. seealso::
      `restart <#restart>`_



restart
+++++++

    - Description: Use restart parameter set to true to start a simulation from scratch (with existing input PDBs and configuration files)
    - Type: ``Boolean``
    - Default: ``False``

    .. seealso::
      `debug <#debug>`_


seed
++++

    - Description: Seed for Platform's pseudo-random numbers generator for reproducibility.
    - Type: ``Integer``
    - Default: ``12345``



