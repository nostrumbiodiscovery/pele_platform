.. _tutorial:

========
Tutorial
========

This section will show how to apply the platform to a Nuclear Hormone real
case

PR Tutorial
============


**Command**::

    python MSM_PELE/main.py MSM_PELE/Examples/PR_1A28_xray_-_minimized.pdb STR Z --cpus 120

**Adaptive Exit Analysis**

An Adaptive Exit simulation is performed and all conformations are
clusterize. In order to check good discretization along the exit
path next instructions should be followed::

	  maestro Pele_STR/PR_1A28_xray_-_minimized_processed.pdb Pele_STR/box.pdb Pele_STR/ouput_clustering/clustering_40_Kmeans.pdb

The representative pdb of each cluster can be found under PELE_STR/output_clustering::

    maestro PELE_STR/ouput_clustering/initial_*

**PELE exploration Analysis**

An exhaustive PELE exploration is performed starting from each cluster from the previous
exit Adaptive period. There are several features we can check to know whether
the exploration was good enough. Run the next commands under Pele_STR/output_pele/

`Transition Study`::

  MSM_PELE/Pele_scripts/AnalysisTools/adaptivePlot 5 6 1000 -rmsd | gnuplot -p

`Total Energy Conformations`::

  MSM_PELE/Pele_scripts/AnalysisTools/adaptivePlot 2 4 1000 -zcol 2 -be
  | guplot -p

`Counter Plot`::

  MSM_PELE/Pele_scripts/AnalysisTools/counter.py 6 10

`Best Binding Energy Structures`::

  MSM_PELE/Pele_scripts/AnalysisTools/bestStructs.py Binding Energy

**MSM Analysis**

`Correlation Analysis`::

  python ~/repos/PyTools/PyTools/autoCorr.py -l 1000 -n 200 --clusters
  MSM_0/clusterCenters_0.dat --trajs allTrajs/

`Correlation Clusters`::

  vmd Pele_STR/output_pele/representative_structures/clusters_0.pdb (Colored
  with the beta property from red (high correlation) to blue (low correlation)

