API
===

Analysis
---------

We included an option to run analysis as API for those of our users, who are familiar with Python. All you have to do is
initialize the Analysis class with the three mandatory parameters (resname, chain and simulation_output) and any optional
you might want to include.

Documentation
+++++++++++++

.. autoclass:: pele_platform.analysis.Analysis
  :members:

Example
++++++++

To start off, you need to initialize Analysis, providing at least the following three parameters:

    * residue name
    * chain ID
    * path to the output folder.

.. code-block:: python

 >> from pele_platform.analysis import Analysis
 >> analysis = Analysis(resname="LIG", chain="Z", simulation_output="LIG_Pele/output")

Once it has been generated, you can call any of the available methods to get the pieces of analysis that interest you or simply
run the whole workflow, which includes extraction of top poses, plots, clustering and PDF report:

.. code-block:: python

    >> analysis.generate(path="my_folder", clustering_type="gaussianmixture", analysis_nclust=3)
