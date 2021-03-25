Analysis
============

Output files
-----------------

Folder structure
++++++++++++++++++

Each simulation will create a number of output folders, as indicated on the tree below.

.. code:: console

    .
    └── LIG_Pele
        ├── output
        └── results
            ├── top_poses
            ├── plots
            └── clusters


The ``output`` folder contains raw output files such as detailed metrics reports and trajectory snapshots for every step, whereas the
``results`` directory holds a more user-friendly, curated output comprising of three separate folders:

- ``top_poses`` - top 100 lowest binding energy structures
- ``plots`` - plots of multiple metrics selected by the user
- ``clusters`` - lowest binding energy cluster representatives and clustering plots.

Detailed metrics
++++++++++++++++++

Additionally, the simulation will create a two CSV files with more detailed metrics:

- ``results/data.csv`` contains a summary of all created poses together with their metrics and the clusters they belong to
- ``results/clusters/info.csv`` provides detailed metric on each cluster, such as its population, mean RMSD, energy percentiles, etc.


Plots
---------
The software will automatically create scatter plots for all metrics, however, if you want to enhance your analysis, you
can check out our interactive and KDE plots.

Interactive plots
+++++++++++++++++++
You can also create your own interactive plots. Simply go to the ``output`` folder and run the following command:

.. code-block:: console

 python -m pele_platform.analysis.interactive_plot 6 5

The script requires two integer arguments indicating the numbers of report columns you would like to plot, in this
case we used columns 5 and 6 corresponding to the binding energy and SASA of the ligand. You can click on data points to get the file names of the structures.

.. image:: ../img/interactive_plot.png
  :width: 400
  :align: center

For more advanced interactive plots, please refer to `PELE++ documentation <https://nostrumbiodiscovery.github.io/pele_docs/intro/GeneralAnalysis/GeneralAnalysis.html>`_.

Kernel density estimate plot
++++++++++++++++++++++++++++++

The KDE plots can aid the simulation analysis by visualising the distribution of ligand poses (similarly to a histogram)
in respect to plotted metrics, such as distance between two atoms (atom_dist) or solvent exposed surface area (SASA).
All you have to do is include the ``kde: true`` flag in your input.yaml. Additionally, you can influence the number of
poses plotted using the ``kde_structs`` flag, where the default included 1000 best energy poses.

.. code-block:: yaml

    kde: true
    kde_structs: 200

Example plot of binding energy vs ligand SASA with a KDE including 1000 best energy poses.

.. image:: ../img/kde.png
  :width: 400
  :align: center

Clusters
-----------

The user has a choice between three clustering methods as well as some control over: mean shift, HDBSCAN and Gaussian mixture model.


API
-----

We included an option to run analysis as API for those of our users, who are familiar with Python. All you have to do is
initialize the Analysis class with the three mandatory parameters (resname, chain and simulation_output) and any optional
you might want to include.

class Analysis
++++++++++++++++
        resname : str (mandatory)
            Residue name of the ligand, e.g. "LIG"
        chain : str (mandatory)
            Chain ID of the ligand, e.g. "Z."
        simulation_output : str (mandatory)
            Path to the output folder of the simulation, e.g. "LIG_Pele/output"
        be_column : int (optional, default = 4)
            Report column with energy metric.
        limit_column : int (optional, default = None)
            Integer specifying the first column from which the meaningful metrics start, e.g. SASA or RMSD.
        traj : str (optional, default = "tarjectory.pdb")
            Trajectory name defaults to "trajectory.pdb", but you should use "trajectory.xtc" if using XTC format.
        report : str (optional, default = "report")
            Report file name, if not using default.
        skip_initial_structures : bool (optional, default = False)
            Skips initial structures (step 0 of the simulation). Should be set to False when running test
            with only one step.
        kde : bool (optional, default = False)
            Set to True to create kernel density estimator plots.
        kde_structs : int (optional, default = 1000)
            Maximum number of structures to consider for the KDE plot.
        topology : str (optional, default = None)
            Path to the topology file, if using XTC trajectories.
        cpus: int (optional, default = 1)
            Number of CPUs to use.

.. code-block:: python

     >> analysis = Analysis(
            resname="LIG",
            chain="Z",
            simulation_output="LIG_Pele/output",
        )

Then you can use one of the available methods to generate top poses, perform clustering or run the whole analysis workflow, e.g.

generate()
++++++++++++

.. code-block:: python

    >> analysis.generate(working_folder, "gaussianmixture")


generate_clusters()
++++++++++++++++++++
Runs the full analysis workflow (plots, top poses and clusters) and saves the results in the supplied path.

        path : str
            The path where the analysis results will be saved
        clustering_type : str (optional, default = 'meanshift')
            The clustering method that will be used to generate the clusters. One of ['gaussianmixture', 'meanshift', 'hdbscan'].
        bandwidth : float (optional, default = 2.5)
            Bandwidth for the mean shift and HDBSCAN clustering.
        analysis_nclust : int (optional, default = 10)
            Number of clusters to create when using the Gaussian mixture model.
        max_top_clusters : int (optional, default = 8)
            Maximum number of clusters to return.
        min_population : float (optional, default = 0.01)
            The minimum amount of structures in a cluster, takes a value between 0 and 1, where 0.01 refers to 1% of all structures.

.. code-block:: python

    >> analysis.generate_clusters(path="my_folder", clustering_type="gaussianmixture", analysis_nclust=3)

generate_plots()
++++++++++++++++++

.. code-block:: python

    >> analysis.generate_plots(working_folder, "gaussianmixture")

generate_report()
++++++++++++++++++

.. code-block:: python

    >> analysis.generate_report(working_folder, "gaussianmixture")

generate_top_poses()
++++++++++++++++++

.. code-block:: python

    >> analysis.generate_top_poses(working_folder, "gaussianmixture")
