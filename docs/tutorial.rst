.. _tutorial:

========
Tutorial
========

Command line tutorial PRR:
--------------------------

python MSM_PELE/main.py MSM_PELE/Examples/PR_1A28_xray_-_minimized.pdb STR Z --cpus 60

Slurm tutorial PRR:
-------------------

- Inside the slurm file change the fields:

	<PDB>: PR_1A28_minimized.pdb
    <resname>: STR
    <chain>: Z

- Launch the job with the next command:

	sbatch MSM_PELE/Examples/conf.sl
