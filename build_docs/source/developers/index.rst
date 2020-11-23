#####################################
Developers: Add your own mode
#####################################


To add your own new simulation mode:


1) Generate your own mode
++++++++++++++++++++++++++++

Add the default of your features below pele_platform/features.py line 27:
  
  
Example:

.. code-block:: python

 "induced_fit_exhaustive" : {"spawning_type": "independent", "bias_column": 5, "epsilon":0.25, "density": "null",
          "simulation_type": "pele", "iterations": 1, "pele_steps": 1000,
          "cluster_values": "[1.75, 2.5, 4, 6]", "cluster_conditions": "[1, 0.6, 0.4, 0.0]",
          "steric_trials": 500, "overlap_factor": 0.65, "params": pcs.INDUCED_FIT,
          "box_radius": 6},

2) Include your mode into the pipeline
++++++++++++++++++++++++++++++++++++++++

Set the type of the simulation below pele_platform/features.py line 90:

Example:

.. code-block:: python

    elif args.induced_fit_fast:
        type_simulation = "induced_fit_fast"

3) Create your mode's input option
+++++++++++++++++++++++++++++++++++

Under pele_platform/main.py line 204 create your input.yaml option. 
Set the default to false so it will get activated when the user specifies induced_fit_fast: true.

Example:

.. code-block:: python

  self.induced_fit_fast = data.get("induced_fit_fast", False)


4) Run your mode
+++++++++++++++++++++

To run your simulation mode in the input file specify the previously generated yaml option.

Example:

.. code-block:: yaml

  induced_fit_fast: true
