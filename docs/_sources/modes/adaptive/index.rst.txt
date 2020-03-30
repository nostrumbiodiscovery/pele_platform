PELE Modes
######################

Automatically configures all simulation parameters to match the desired configuration. **To choose between:
induce fit, pocket exploration, binding path exploration, minimization, bias exploration, exit path**


Induced fit_fast
==================

- **induced_fit_fast**: Run a short induced fit AdaptivePELE simulation paramaters by setting the center of the box in the
  cm of the ligand, a box radius of 10A, small rotations and translations and a high number of
  steric clashes and sidechain predition frequency. **Useful to refine docking poses, and search
  new conformations within the same binding site**.

..  code-block:: yaml

  induced_fit_fast: true

Induced fit_exhaustive
========================

- **induced_fit_exhaustive**: Run a long induced fit PELE simulation paramaters by setting the center of the box in the
  cm of the ligand, a box radius of 10A, small rotations and translations and a high number of
  steric clashes and sidechain predition frequency. **Useful to refine docking poses, and search
  new conformations within the same binding site**.

..  code-block:: yaml

  induced_fit_exhaustive: true


Pocket Exploration
=====================

- **ppi**: Configure a global exploration by randomizing the ligand all around the protein. Then the simulation will start from all configurationsof the system at the same time. The number of configurations (ligand-protein systems) can be chosen thorugh the **poses** flag by default will be cpus-1. Once finished, an exhaustive induced fit exploration will follow to retrieve all best identified pockets and putative binding modes. **Useful to search for putative binding modes when neither pose or pocket are known.**

..  code-block:: yaml

  ppi: true


Rescoring
============

Simulation to refine around an initial conformation. Not looking to find a new binding mode but to minimize
the actual one. **Useful to minimize certain pose or to perform data augmentation for ML/DL methods**

..  code-block:: yaml

  rescoring: true

Local Exploration
=====================

- **out_in**: Local exploration to move the ligand from the bulk to the binding site. The box center set on the
  center of mass of the ligand with a radius of 30A, steering 1 50% of the times, and a slight bias towards binding energies.
  **Useful when no docking is possible in the binding site and you need to open up the pocket.**

..  code-block:: yaml

  out_in: true

Biased
=========

- **bias**: Bias exploration towards the indicated bias column. The box center is set on the center of mass of the ligand with
  a radius of 30A, and a bias towards the chosen metric is set. An epsilon fraction of processors are distributed proportionally to the value of a metric, and the rest are inverselyProportional distributed. Therefore, the **epsilon** value controls fraction of the processors that will be assigned according to the selected metric in **biascolumn**. **Useful to speed up local explorations when the desired binding site is known**


..  code-block:: yaml

  epsilon: 0.5
  bias_column: 5 (starting by 1 on the reports)

Exit path
==============

- **in_out**: Explore the dissociative path of a molecule. At each step the box is center on the most exterior cluster
  and there is a bias towards higher values of SASA. This type accepts a **exit_metric** which represents a column in the report file, an **exit_value** which represents a value for the metric and a **exit_condition** parameter which can be either “<” or “>”, default value is “<”. The simulation will terminate after the metric written in the metricCol reaches a value smaller or greater than exitValue, depending on the condition specified. An example of the exit condition block that would terminate the program after 4 trajectories reaches a value of more than 0.9 for the sixth column (6th starting to count from 1) of the report file would look like below. **Useful when studying the dissoative binding path, or to open up the pocket**


..  code-block:: yaml

  in_out: true
  exit_value: 0.9
  exit_condition: ">"
  exit_trajnum: 4

