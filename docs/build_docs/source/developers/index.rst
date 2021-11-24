========================
Advanced tips and tricks
========================


Enable PELE InterStepLogger
---------------------------------
InterStepLogger functionality allows capturing structural and simulation data during different stages of a step.
It is useful for debugging purposes.

To enable this option you need to set **inter_step_logger** parameter to true in your **input.yaml**:

..  code-block:: yaml

  inter_step_logger: true

The new report and trajectories files were saved in the subdirectory ``InterStepLogs``, inside the PELE++ output folder.


Add custom simulation mode
-----------------------------

If your simulation set up requires tweaking a lot of parameters, you can add a new **feature** to the platform.

1. Generate your own mode
^^^^^^^^^^^^^^^^^^^^^^^^^

Add a new feature with your own defaults under ``simulation_params`` key in ``pele_platform/features/adaptive.py``, e.g.

.. code-block:: python

 "my_custom_induced_fit" : {"spawning_type": "independent", "bias_column": 5, "epsilon":0.25, "density": "null",
          "simulation_type": "pele", "iterations": 5, "pele_steps": 500,
          "cluster_values": "[1.75, 2.5, 4, 6]", "cluster_conditions": "[1, 0.6, 0.4, 0.0]",
          "steric_trials": 500, "overlap_factor": 0.65, "params": pcs.INDUCED_FIT,
          "box_radius": 6},

2. Include it in the pipeline
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Assign the type of the simulation depending on the YAML arguments in ``pele_platform/features/adaptive.py``, e.g.

.. code-block:: python

    elif args.my_custom_induced_fit:
        type_simulation = "my_custom_induced_fit"

3. Add an input option
^^^^^^^^^^^^^^^^^^^^^^

A new input flag will not be recognised by the platform, unless you add it to:

    - ``_parse`` function in ``pele_platform/Utilities/Helpers/yaml_parser.py`` and set the default to None
    - ``VALID_FLAGS_PLATFORM`` dictionary in ``pele_platform/Checker/valid_flags.py``

Example in ``yaml_parser.py``:

.. code-block:: python

  self.induced_fit_fast = data.get(valid_flags["my_custom_induced_fit"], None)

4. Run your mode
^^^^^^^^^^^^^^^^

To run your simulation mode in the input file, include the previously generated YAML flag, for example:

.. code-block:: yaml

  induced_fit_fast: true
