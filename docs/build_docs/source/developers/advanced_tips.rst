Advanced tips & tricks
============================

1. Enable PELE InterStepLogger
---------------------------------
InterStepLogger functionality allows capturing structural and simulation data during different stages of a step.
It is useful for debugging purposes.

To enable this option you need to set **inter_step_logger** parameter to true in your **input.yaml**:

..  code-block:: yaml

  inter_step_logger: true

The new report and trajectories files were saved in the subdirectory ``InterStepLogs``, inside the PELE++ output folder. 
