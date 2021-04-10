"""
This module manages the analysis toolkit of the platform.
"""

__all__ = ["Analysis"]


class Analysis(object):
    """
    General class to manage all analysis operations.
    """

    _EPOCH_LABEL = "epoch"
    _TRAJECTORY_LABEL = "trajectory"
    _REPORT = "report"
    _STEP_LABEL = "numberOfAcceptedPeleSteps"

    def __init__(self, resname, chain, simulation_output,
                 be_column=4, limit_column=None, traj="trajectory.pdb",
                 report=None, skip_initial_structures=True, kde=False,
                 kde_structs=1000, topology=None, cpus=1):
        """
        It initializes an Analysis instance which it depends on
        the general Parameters class of the PELE Platform.

        Parameters
        ----------
        resname : str
            Residue name of the ligand, e.g. "LIG"
        chain : str
            Chain ID of the ligand, e.g. "Z."
        simulation_output : str
            Path to the output folder of the simulation, e.g.
            "LIG_Pele/output"
        be_column : int
            Column with energy metric, default 4.
        limit_column : int
            Integer specifying the first column from which the meaningful
            metrics start, e.g. SASA or RMSD.
        traj : str
            Trajectory name defaults to "trajectory.pdb",
            but you should use "trajectory.xtc" if using XTC format.
        report : str
            Report file name, if not using default.
        skip_initial_structures : bool
            Skips initial structures (step 0 of the simulation),
            default is True. Should be set to True when running test
            with only one step
        kde : bool
            Set to True to create kernel density estimator plots. Default
            is False
        kde_structs : int
            Maximum number of structures to consider for the KDE plot.
            Default is 1000
        topology : str
            Path to the topology file, if using XTC trajectories. Default
             is None
        cpus : int
            Number of CPUs to use. Default is 1
        """
        from pele_platform.analysis import DataHandler

        self.residue = resname
        self.chain = chain
        self.output = simulation_output
        self.be_column = be_column
        if be_column is None:
            self.be_column = 4
        self.limit_column = limit_column
        self.kde = kde
        self.kde_structs = kde_structs
        self.traj = traj
        self.report = report if report else self._REPORT
        self.skip_initial_structures = skip_initial_structures
        self.topology = topology
        self.cpus = cpus

        self._data_handler = DataHandler(
            sim_path=self.output,
            report_name=self.report,
            trajectory_name=self.traj,
            be_column=self.be_column,
            skip_initial_structures=self.skip_initial_structures)
        self._dataframe = self._data_handler.get_reports_dataframe()

    @classmethod
    def from_parameters(cls, parameters):
        """
        It initializes an Analysis object from a Parameters object.

        Parameters
        ----------
        parameters : a Parameters object
            The Parameters object containing the parameters that belong
            to the simulation

        Returns
        -------
        analysis : an Analysis object
            The Analysis object obtained from the parameters that were
            supplied
        """
        import os

        simulation_output = os.path.join(parameters.pele_dir, parameters.output)

        analysis = Analysis(resname=parameters.residue,
                            chain=parameters.chain,
                            simulation_output=simulation_output,
                            be_column=parameters.be_column,
                            limit_column=parameters.limit_column,
                            traj=parameters.traj_name,
                            report=parameters.report_name,
                            skip_initial_structures=not parameters.test,
                            kde=parameters.kde,
                            kde_structs=parameters.kde_structs,
                            topology=parameters.topology,
                            cpus=parameters.cpus)

        return analysis

    @property
    def parameters(self):
        """
        It returns the attributes of this Analysis object as a dictionary.

        Returns
        -------
        params : dict
            A dictionary of parameters
        """
        params = {key: value for key, value in self.__dict__.items()
                  if key[:1] != "_"}
        return params

    def get_dataframe(self, filter=False, threshold=None):
        """
        Parameters
        ----------
        filter : bool
            Whether to filter the entries with highest energies according
            to the threshold value or not. Default is False
        threshold : float
            The ratio of high-energy entries that will be filtered out.
            Default is None and will be initialized with a threshold of
            0.02

        Returns
        -------
        dataframe : a pandas.DataFrame object
            The dataframe containing the information from PELE reports
        """
        if filter:
            return self._data_handler.remove_outliers_from_dataframe(
                self._dataframe, threshold)
        else:
            return self._dataframe

    def dataframe_to_csv(self, path, filter=False, threshold=None):
        """
        It saves the dataframe in the supplied path as a csv file.

        Parameters
        ----------
        path : str
            The path where the dataframe will be saved
        filter : bool
            Whether to filter the entries with highest energies according
            to the threshold value or not. Default is False
        threshold : float
            The ratio of high-energy entries that will be filtered out.
            Default is None and will be initialized with a threshold of
            0.02
        """
        # Get dataframe
        dataframe = self.get_dataframe(filter, threshold)

        # Save it as a csv file
        dataframe.to_csv(path, index=False)

    def generate(self, path, clustering_type='meanshift',
                 bandwidth=2.5, analysis_nclust=10,
                 max_top_clusters=8,
                 top_clusters_criterion="interaction_25_percentile",
                 min_population=0.01, max_top_poses=100,
                 representatives_criterion="interaction_5_percentile"):
        """
        It runs the full analysis workflow (plots, top poses and clusters)
        and saves the results in the supplied path.

        Parameters
        ----------
        path : str
            The path where the analysis results will be saved
        clustering_type : str
            The clustering method that will be used to generate the
            clusters. One of ['gaussianmixture', 'meanshift', 'hdbscan'].
            Default is 'meanshift'
        bandwidth : float
            Bandwidth for the mean shift and HDBSCAN clustering. Default is
            2.5
        analysis_nclust : int
            Number of clusters to create when using the Gaussian mixture
            model. Default is 10
        max_top_clusters : int
            Maximum number of clusters to return. Default is 8
        min_population : float
            The minimum amount of structures in a cluster, takes a value
            between 0 and 1. Default is 0.01 (i.e. 1%)
        max_top_poses : int
            Number of top poses to retrieve. Default = 100.
        top_clusters_criterion : str
            Criterion to select top clusters. Default is
            "interaction_25_percentile". One of ["total_25_percentile",
            "total_5_percentile", "total_mean", "total_min",
            "interaction_25_percentile", "interaction_5_percentile",
            "interaction_mean", "interaction_min", "population"]
        representatives_criterion : str
            Criterion to select cluster representative structures. Default is
            "interaction_5_percentile". One of ["total_25_percentile",
            "total_5_percentile", "total_mean", "total_min",
            "interaction_25_percentile", "interaction_5_percentile",
            "interaction_mean", "interaction_min"]
        """
        import os
        path = self._check_existing_directory(path)

        summary_file = os.path.join(path, "data.csv")
        plots_folder = os.path.join(path, "plots")
        top_poses_folder = os.path.join(path, "top_poses")
        clusters_folder = os.path.join(path, "clusters")
        report_file = os.path.join(path, "summary.pdf")

        if not os.path.exists(plots_folder):
            os.makedirs(plots_folder)
        if not os.path.exists(top_poses_folder):
            os.makedirs(top_poses_folder)
        if not os.path.exists(clusters_folder):
            os.makedirs(clusters_folder)

        # Save dataframe
        # (it will be save later again, replacing this file, to include
        # the column of clusters, only if clustering can run successfully)
        self._dataframe.to_csv(summary_file, index=False)

        # Generate analysis results
        self.generate_plots(plots_folder)
        best_metrics = self.generate_top_poses(top_poses_folder, max_top_poses)
        self.generate_clusters(clusters_folder, clustering_type,
                               bandwidth, analysis_nclust,
                               max_top_clusters, top_clusters_criterion,
                               min_population, representatives_criterion)

        self.generate_report(plots_folder, top_poses_folder,
                             clusters_folder, best_metrics, report_file)

    def generate_plots(self, path):
        """
        It generates the plots.

        Parameters
        ----------
        path : str
            The path where the plots will be saved
        """
        from pele_platform.analysis import Plotter

        # Get dataframe, filtering highest 2% energies out
        dataframe = self.get_dataframe(filter=True)

        # Initialize plotter
        plotter = Plotter(dataframe)

        metrics = self._data_handler.get_metrics()

        # In case there is a Interaction energy we will generate two different
        # plots for each metric we find:
        #  - Total energy vs Interaction energy vs Metric (in the color bar)
        #  - Interaction energy vs Metric
        # In case there Interaction energy is not available, we will generate
        # a single plot for each metric we find:
        #  - Total energy vs Metric
        if "Binding Energy" in metrics:
            t_energy = "currentEnergy"
            i_energy = "Binding Energy"
        else:
            t_energy = "currentEnergy"
            i_energy = None

        # The minimum value for the limit column is 4, since previous
        # columns in PELE report files does not contain any metric
        if self.limit_column is not None and self.limit_column > 4:
            limit_column = self.limit_column - 4
        else:
            limit_column = 0

        # Iterate over all the metrics found in the reports
        for metric in metrics[limit_column:]:
            # Avoid comparing an energy with itself
            if metric == t_energy or metric == i_energy:
                continue

            if i_energy is not None:
                plotter.plot_two_metrics(t_energy, i_energy, metric,
                                         output_folder=path)
                plotter.plot_two_metrics(metric, i_energy,
                                         output_folder=path)
            else:
                plotter.plot_two_metrics(metric, t_energy,
                                         output_folder=path)

            if self.kde:
                plotter.plot_kde(metric, i_energy, output_folder=path,
                                 kde_structs=self.kde_structs)

    def generate_top_poses(self, path, n_poses):
        """
        It selects and saves the top poses.

        Parameters
        ----------
        path : str
            The path where the top poses will be saved
        n_poses : int
            The number of top poses to retrieve

        Returns
        -------
        best_metrics : list[float]
            The list that contains the metrics belonging to the extracted
            best poses
        """
        # Get metrics and locate Interaction energy
        metrics = self._data_handler.get_metrics()

        if "Binding Energy" in metrics:
            metric = "Binding Energy"
        else:
            metric = "currentEnergy"

        print("Retrieve {} Best Poses".format(n_poses))

        top_poses = self._data_handler.get_top_entries(metric, n_poses)
        best_metrics = self._extract_poses(top_poses, metric, path)

        return best_metrics

    def generate_clusters(self, path, clustering_type,
                          bandwidth=2.5, analysis_nclust=10,
                          max_top_clusters=8,
                          top_clusters_criterion="interaction_25_percentile",
                          min_population=0.01,
                          representatives_criterion="interaction_5_percentile"):
        """
        It generates the structural clustering of ligand poses.

        Parameters
        ----------
        path : str
            The path where the clusters will be saved
        clustering_type : str
            The clustering method that will be used to generate the
            clusters
        bandwidth : float
            Bandwidth for the mean shift and HDBSCAN clustering. Default is
            2.5
        analysis_nclust : int
            Number of clusters to create when using the Gaussian mixture
            model. Default is 10
        max_top_clusters : int
            Maximum number of clusters to return. Default is 8
        top_clusters_criterion : str
            Criterion to select top clusters. Default is
            "interaction_25_percentile". One of ["total_25_percentile",
            "total_5_percentile", "total_mean", "total_min",
            "interaction_25_percentile", "interaction_5_percentile",
            "interaction_mean", "population"]
        min_population : float
            The minimum amount of structures in a cluster, takes a value
            between 0 and 1. Default is 0.01 (i.e. 1%)
        representatives_criterion : str
            Criterion to select cluster representative structures. Default is
            "interaction_5_percentile". One of ["total_25_percentile",
            "total_5_percentile", "total_mean", "total_min",
            "interaction_25_percentile", "interaction_5_percentile",
            "interaction_mean", "interaction_min"]
        """
        import os
        from pele_platform.Utilities.Helpers.helpers import check_make_folder
        from pele_platform.constants.constants import \
            metric_top_clusters_criterion, cluster_representatives_criterion

        check_make_folder(path)

        # Get clustering object
        clustering, max_coordinates = self._get_clustering(clustering_type,
                                                           bandwidth,
                                                           analysis_nclust)

        # Extract coordinates
        coordinates, dataframe = self._extract_coordinates(max_coordinates)

        # Skip clustering in case
        if coordinates is None or dataframe is None:
            return

        # Filter coordinates
        coordinates, dataframe, energetic_threshold = \
            self._filter_coordinates(coordinates, dataframe)

        # Cluster coordinates
        print(f"Cluster ligand binding modes")
        clusters = clustering.get_clusters(coordinates, self._dataframe,
                                           dataframe, os.path.dirname(path))

        # Analyze and save clustering results
        rmsd_per_cluster = self._calculate_cluster_rmsds(clusters, coordinates)
        cluster_summary = self._analyze_clusters(clusters, dataframe,
                                                 rmsd_per_cluster)

        if len(cluster_summary) == 0:
            print(f"No clusters could be obtained, " +
                  f"clustering analysis is skipped")

            return

        cluster_subset, cluster_summary = \
            self._select_top_clusters(clusters, cluster_summary,
                                      top_clusters_criterion,
                                      max_clusters_to_select=max_top_clusters,
                                      min_population_to_select=min_population)
        print(f"Retrieve top clusters based on " +
              f"{metric_top_clusters_criterion[top_clusters_criterion]}.")

        # Save cluster summary to file with information about selected labels
        cluster_summary.to_csv(os.path.join(path, "info.csv"), index=False)

        # Plot cluster descriptors
        self._plot_cluster_descriptors(cluster_subset, dataframe,
                                       cluster_summary, path)

        # Plot clusters
        self._plot_clusters(cluster_subset, dataframe, cluster_summary, path)

        # Save cluster representative structures
        self._save_cluster_representatives(cluster_subset, dataframe, path,
                                           representatives_criterion)
        print(
            f"Retrieve top cluster representative structures based on " +
            f"{cluster_representatives_criterion[representatives_criterion]}.")

    def generate_report(self, plots_path, poses_path, clusters_path,
                        best_metrics, filename):
        """
        It generates the final simulation report as a PDF file.

        Parameters
        ----------
        plots_path : str
            The path where the plots are saved
        poses_path : str
            The path where the top poses are saved
        clusters_path : str
            The path where the clusters are saved
        best_metrics : list[float]
            The list that contains the metrics belonging to the extracted
            best poses
        filename : str
            The filename for the simulation report
        """
        import os
        import glob
        from pele_platform.analysis import pdf_report

        plots = glob.glob(os.path.join(plots_path, "*.png"))
        poses = glob.glob(os.path.join(poses_path, "*"))
        clusters = glob.glob(os.path.join(clusters_path, "*.png"))

        report = pdf_report.create_report(plots, clusters, poses,
                                          best_metrics, filename)

        print("PDF summary report successfully written to: {}".format(report))

    def _get_clustering(self, clustering_type, bandwidth, analysis_nclust):
        """
        It returns the clustering object according to the supplied
        clustering type.

        Parameters
        ----------
        clustering_type : str
            The type of clustering to use
        bandwidth : float
            Bandwidth for the mean shift and HDBSCAN clustering. Default is
            2.5
        analysis_nclust : int
            Number of clusters to create when using the Gaussian mixture
            model. Default is 10

        Returns
        -------
        clustering : a Clustering object
            The Clustering object that matches with the supplied
            clustering type
        max_coordinates : int
            The maximum number of coordinates to extract from the
            residue
        """
        from pele_platform.analysis import (GaussianMixtureClustering,
                                            HDBSCANClustering,
                                            MeanShiftClustering)

        if clustering_type.lower() == "gaussianmixture":
            clustering = GaussianMixtureClustering(analysis_nclust)
            max_coordinates = 10
        elif clustering_type.lower() == "hdbscan":
            clustering = HDBSCANClustering(bandwidth)
            max_coordinates = 10
        elif clustering_type.lower() == "meanshift":
            clustering = MeanShiftClustering(bandwidth)
            max_coordinates = 5
        else:
            raise ValueError("Invalid clustering type: " +
                             "'{}'. ".format(clustering_type) +
                             "It should be one of ['GaussianMixture', " +
                             "'HDBSCAN', 'MeanShift']")

        return clustering, max_coordinates

    def _extract_coordinates(self, max_coordinates):
        """
        It extracts the coordinates of the simulation and creates the
        dataframe with the metrics of each snapshot.

        Parameters
        ----------
        max_coordinates : int
            The maximum number of coordinates to extract from the
            residue

        Returns
        -------
        coordinates : numpy.array
            The array of coordinates to filter
        dataframe : a pandas.dataframe object
            The dataframe containing the PELE reports information that
            follows the same ordering as the array of coordinates
        """
        print(f"Extract coordinates for clustering")

        if not self.topology:
            coordinates, dataframe = self._data_handler.extract_PDB_coords(
                self.residue, remove_hydrogen=True,
                n_proc=self.cpus, max_coordinates=max_coordinates)
        else:
            coordinates, dataframe = self._data_handler.extract_XTC_coords(
                self.residue, self.topology,
                remove_hydrogen=True, max_coordinates=max_coordinates)

        if coordinates is None or dataframe is None:
            print(f"Coordinate extraction failed, " +
                  f"clustering analysis is skipped")
            return None, None

        if len(coordinates) < 2:
            print(f"Not enough coordinates, " +
                  f"clustering analysis is skipped")
            return None, None

        return coordinates, dataframe

    def _extract_poses(self, dataframe, metric, output_path):
        """
        Given a dataframe, it extracts all the corresponding poses.

        Parameters
        ----------
        dataframe : a pandas.DataFrame object
            The dataframe containing the information of the poses to extract
        metric : str
            The metric to highlight in the filename
        output_path : str
            The path where the poses will be extracted

        Returns
        -------
        values : list[float]
            The list that contains the metrics belonging to the extracted
            best poses
        """
        import os
        from pele_platform.Utilities.Helpers import get_suffix
        from pele_platform.Utilities.Helpers.bestStructs import (
            extract_snapshot_from_pdb,
            extract_snapshot_from_xtc)

        values = dataframe[metric].tolist()
        paths = dataframe[self._TRAJECTORY_LABEL].tolist()
        epochs = dataframe[self._EPOCH_LABEL].tolist()
        file_ids = [get_suffix(os.path.splitext(trajectory)[0])
                    for trajectory in paths]
        steps = list(map(int, dataframe[self._STEP_LABEL].tolist()))

        # To prevent hiding files in case epochs is a list of empty strings
        if all([len(epoch) == 0 for epoch in epochs]):
            epochs = [0, ] * len(values)

        # TODO which is the purpose of this hardcoded distance_key?
        distance_key = "distance0.5"
        if distance_key in dataframe.columns:
            dist_values = dataframe[distance_key].tolist()
            filename_template = "{}.{}.{}_BindEner{:.2f}_AtomDist{:.2f}.pdb"
            file_names = \
                [filename_template.format(epoch, report, step, value, dist)
                 for epoch, step, report, value, dist
                 in zip(epochs, steps, file_ids, values, dist_values)]
        else:
            filename_template = "{}.{}.{}_BindEner{:.2f}.pdb"
            file_names = [filename_template.format(epoch, report, step, value)
                          for epoch, step, report, value
                          in zip(epochs, steps, file_ids, values)]

        # Read trajectory and output snapshot
        for f_id, f_out, step, path in zip(file_ids, file_names, steps, paths):
            if not self.topology:
                try:
                    extract_snapshot_from_pdb(path=path,
                                              f_id=f_id,
                                              output=output_path,
                                              topology=self.topology,
                                              step=step,
                                              out_freq=1,
                                              f_out=f_out)
                except UnicodeDecodeError:
                    raise Exception("XTC output being treated as PDB. " +
                                    "Please specify XTC with the next " +
                                    "flag. traj: 'trajectory_name.xtc' " +
                                    "in the input.yaml")
            else:
                extract_snapshot_from_xtc(path=path,
                                          f_id=f_id,
                                          output=output_path,
                                          topology=self.topology,
                                          step=step,
                                          out_freq=1,
                                          f_out=f_out)

        return values

    def _analyze_clusters(self, clusters, dataframe, rmsd_per_cluster):
        """
        It analyzes the clusters and generates a summary with all
        the calculated descriptors. It also generates some plots
        with information about the clusters.

        Parameters
        ----------
        clusters : a numpy.array object
            The array of cluster labels that were obtained
        dataframe : a pandas.dataframe object
            The dataframe containing the PELE reports information that
            follows the same ordering as the array of clusters
        rmsd_per_cluster : dict[int, float]
            The mean RMSD of each cluster

        Returns
        -------
        cluster_summary : a pandas.dataframe object
            The dataframe containing summary of all clusters that were
            analyzed
        """
        from collections import defaultdict
        import pandas as pd
        import numpy as np

        metrics = self._data_handler.get_metrics()

        clusters_population = defaultdict(int)
        for cluster in clusters:
            clusters_population[cluster] += 1

        summary = pd.DataFrame([(cluster, population / len(clusters),
                                 rmsd_per_cluster[cluster])
                                for cluster, population
                                in clusters_population.items()
                                if not cluster < 0],
                               columns=["Cluster", "Population", "MeanRMSD"])

        # Generate descriptors and boxplots for all reported metrics
        descriptors = defaultdict(dict)
        for metric in metrics:
            values_per_cluster = defaultdict(list)
            values = list(dataframe[metric])

            if len(clusters) != len(values):
                print("Warning: metric '{}' ".format(metric) +
                      "array has a wrong size. It will be skipped " +
                      "from the clustering analysis. Expected size: " +
                      "{}".format(len(clusters)))
                continue

            # Arrange metrics per cluster
            for cluster, value in zip(clusters, values):
                # Skip outliers
                if cluster < 0:
                    continue
                values_per_cluster[cluster].append(value)

            # Calculate descriptors
            for cluster in values_per_cluster:
                descriptors["{} min".format(metric)][cluster] = \
                    np.min(values_per_cluster[cluster])
                descriptors["{} 5-percentile".format(metric)][cluster] = \
                    np.percentile(values_per_cluster[cluster], 5)
                descriptors["{} 25-percentile".format(metric)][cluster] = \
                    np.percentile(values_per_cluster[cluster], 25)
                descriptors["{} mean".format(metric)][cluster] = \
                    np.mean(values_per_cluster[cluster])
                descriptors["{} 75-percentile".format(metric)][cluster] = \
                    np.percentile(values_per_cluster[cluster], 75)
                descriptors["{} 95-percentile".format(metric)][cluster] = \
                    np.percentile(values_per_cluster[cluster], 95)
                descriptors["{} max".format(metric)][cluster] = \
                    np.max(values_per_cluster[cluster])
                descriptors["{} standard deviation".format(metric)][cluster] = \
                    np.std(values_per_cluster[cluster])

        # Add descriptors to summary dataframe
        for label, values_per_cluster in descriptors.items():
            summary[label] = [values_per_cluster[cluster]
                              for cluster in summary["Cluster"]
                              if not cluster < 0]

        return summary

    def _select_top_clusters(self, clusters, cluster_summary,
                             top_clusters_criterion,
                             max_clusters_to_select,
                             min_population_to_select):
        """
        It selects the top clusters based on a user-defined metric (or 5th percentile Binding Energy as a default).

        Parameters
        ----------
        clusters : a numpy.array object
            The array of cluster labels that were obtained
        cluster_summary : a pandas.dataframe object
            The dataframe containing the summary of all clusters that were
            analyzed
        top_clusters_criterion : str
            Criterion to select top clusters. One of ["total_25_percentile",
            "total_5_percentile", "total_mean", "total_min",
            "interaction_25_percentile", "interaction_5_percentile",
            "interaction_mean", "interaction_min", "population"]
        max_clusters_to_select : int
            The maximum number of clusters to select as top
        min_population_to_select : float
            The minimum population the clusters must have in order to
            be selected

        Returns
        -------
        cluster_subset : a numpy.array object
            The array of cluster after the selection. Those clusters
            that were not selected are now labeled with a -1
        cluster_summary : a pandas.dataframe object
            The dataframe containing summary of all clusters that were
            analyzed. It is updated with an extra column containing
            cluster names of top clusters
        """
        from pele_platform.analysis.clustering import get_cluster_label
        from pele_platform.constants.constants import \
            metric_top_clusters_criterion

        # Get metric to be used in the top cluster selection
        if top_clusters_criterion.lower() in metric_top_clusters_criterion:
            user_metric = \
                metric_top_clusters_criterion[top_clusters_criterion.lower()]
        else:
            raise ValueError('Invalid top_clusters_criterion ' +
                             '\'{}\''.format(top_clusters_criterion.lower()) +
                             '. It must be one of ' +
                             '{}'.format(metric_top_clusters_criterion.keys()))

        # Check if the selected metric is available
        if user_metric in list(cluster_summary.columns):
            metric = user_metric
        else:
            print('Warning: supplied metric for the top cluster selection ' +
                  'is missing in the reports, '
                  '\'{}\'. '.format(top_clusters_criterion) +
                  'Cluster population will be used instead.')
            metric = "Population"

        # Filter cluster summary by Population
        filtered_cluster_summary = \
            cluster_summary[cluster_summary["Population"] >=
                            min_population_to_select]

        if len(filtered_cluster_summary) == 0:
            print('Warning: no cluster fulfills the minimum population '
                  'threshold. Consider increasing the cluster size or ' +
                  'lowering the minimum population value.')

        # Select top clusters based on the chosen metric
        if metric == "Population":
            filtered_cluster_summary = \
                filtered_cluster_summary.nlargest(max_clusters_to_select,
                                                  metric)
        else:
            filtered_cluster_summary = \
                filtered_cluster_summary.nsmallest(max_clusters_to_select,
                                                   metric)

        top_clusters = list(filtered_cluster_summary["Cluster"])

        cluster_reindex_map = {}
        for index, cluster in enumerate(sorted(top_clusters)):
            cluster_reindex_map[cluster] = index

        cluster_subset = []
        for cluster in clusters:
            if cluster in top_clusters:
                cluster_subset.append(cluster_reindex_map[cluster])
            else:
                cluster_subset.append(-1)

        cluster_summary["Selected labels"] = [
            get_cluster_label(cluster_reindex_map[cluster])
            if cluster in cluster_reindex_map
            else "-" for cluster in cluster_summary["Cluster"]]

        return cluster_subset, cluster_summary

    def _plot_cluster_descriptors(self, clusters, dataframe,
                                  cluster_summary, path):
        """
        It plots cluster descriptors.

        Parameters
        ----------
        clusters : a numpy.array object
            The array of cluster labels that were obtained
        dataframe : a pandas.dataframe object
            The dataframe containing the PELE reports information that
            follows the same ordering as the array of clusters
        cluster_summary : a pandas.dataframe object
            The dataframe containing summary of all clusters that were
            analyzed
        path : str
            The path where the output files will be saved at
        """
        import os
        from collections import defaultdict
        from matplotlib import pyplot as plt
        from pele_platform.analysis.clustering import get_cluster_label

        if not os.path.exists(path):
            os.mkdir(path)

        sorted_summary = cluster_summary.sort_values(by=['Cluster'],
                                                     inplace=False,
                                                     ascending=True)

        xticks = list()
        xticklabels = list()
        for cluster_id, cluster_label in zip(sorted_summary['Cluster'],
                                             sorted_summary['Selected labels']):
            if cluster_label != '-':
                xticks.append(cluster_id)
                xticklabels.append(cluster_label)
            elif cluster_id % 10 == 0:
                if len(xticks) > 0:
                    if abs(xticks[-1] - cluster_id) > 5:
                        xticks.append(cluster_id)
                        xticklabels.append(cluster_id)
                else:
                    xticks.append(cluster_id)
                    xticklabels.append(cluster_id)

        # Plot Mean RMSD per cluster
        fig, ax = plt.subplots()
        ax.scatter(sorted_summary['Cluster'],
                   sorted_summary['MeanRMSD'],
                   s=[population * 300 for population
                      in sorted_summary['Population']],
                   label="Cluster population")
        ax.set_xlabel("Cluster label/id")
        ax.set_ylabel("Mean RMSD (Å)")
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticklabels)
        plt.legend()
        plot_filename = os.path.join(path, "clusters_meanRMSD.png")
        plt.savefig(plot_filename)

        metrics = self._data_handler.get_metrics()
        for metric in metrics:
            values_per_cluster = defaultdict(list)
            values = list(dataframe[metric])

            if len(clusters) != len(values):
                print("Warning: metric '{}' ".format(metric) +
                      "array has a wrong size. It will be skipped " +
                      "from the clustering analysis. Expected size: " +
                      "{}".format(len(clusters)))
                continue

            # Arrange metrics per cluster
            for cluster, value in zip(clusters, values):
                # Skip outliers
                if cluster < 0:
                    continue
                values_per_cluster[cluster].append(value)

            # Generate boxplots
            try:
                fig, ax = plt.subplots()

                ax.boxplot([values_per_cluster[cluster]
                            for cluster in sorted(values_per_cluster.keys())])
                ax.set_xticklabels([get_cluster_label(cluster)
                                    for cluster
                                    in sorted(values_per_cluster.keys())])

                ax.set_ylabel(metric)
                ax.set_xlabel("Cluster label")

                boxplot_filename = \
                    os.path.join(path,
                                 "top_clusters_{}_boxplot.png".format(metric))
                boxplot_filename.replace(' ', '_')

                plt.savefig(boxplot_filename)

            except IndexError:
                print("Samples too disperse to produce a cluster " +
                      "for metric {}".format(metric))

    def _plot_clusters(self, clusters, dataframe, cluster_summary, path):
        """
        It generates the plots for clusters.

        Parameters
        ----------
        clusters : a numpy.array object
            The array of cluster labels that were obtained
        dataframe : a pandas.dataframe object
            The dataframe containing the PELE reports information that
            follows the same ordering as the array of clusters
        cluster_summary : a pandas.dataframe object
            The dataframe containing summary of all clusters that were
            analyzed
        path : str
            The path where the output file will be saved at
        """
        from pele_platform.analysis import Plotter

        metrics = self._data_handler.get_metrics()

        # Initialize plotter
        plotter = Plotter(dataframe)

        if "Binding Energy" in metrics:
            energy = "Binding Energy"
            plotter.plot_clusters("currentEnergy", energy,
                                  output_folder=path, clusters=clusters)
        else:
            energy = "currentEnergy"

        # The minimum value for the limit column is 4, since previous
        # columns in PELE report files does not contain any metric
        if self.limit_column is not None and self.limit_column > 4:
            limit_column = self.limit_column - 4
        else:
            limit_column = 0

        # Iterate over all the metrics found in the reports
        for metric in metrics[limit_column:]:
            plotter.plot_clusters(metric, energy, output_folder=path,
                                  clusters=clusters)

    def _filter_coordinates(self, coordinates, dataframe, threshold=0.25):
        """
        It filters the coordinates by total energy according to the
        threshold that is supplied. A threshold of 0.25 means that the
        25% of structures with highest energies will be discarded.

        Parameters
        ----------
        coordinates : numpy.array
            The array of coordinates to filter
        dataframe : a pandas.dataframe object
            The dataframe containing the PELE reports information that
            follows the same ordering as the array of coordinates
        threshold : float
            A value between 0 and 1 that defines the ratio of structures
            to discard. Default is 0.25

        Returns
        -------
        filtered_coordinates : numpy.array
            The array of coordinates resulting from the filtering
        filtered_dataframe : a pandas.dataframe object
            The dataframe resulting from the filtering
        energetic_threshold : float
            The energetic value that fulfills the supplied threshold
        """
        import numpy as np

        total_energies = list(dataframe['currentEnergy'])
        energetic_threshold = np.quantile(total_energies, 1 - threshold)

        filtered_coordinates = []
        for coors_array, total_energy in zip(coordinates, total_energies):
            if total_energy <= energetic_threshold:
                filtered_coordinates.append(coors_array)
        filtered_coordinates = np.array(filtered_coordinates)

        filtered_dataframe = \
            dataframe.query('currentEnergy<={}'.format(energetic_threshold))

        return filtered_coordinates, filtered_dataframe, energetic_threshold

    def _calculate_cluster_rmsds(self, clusters, coordinates):
        """
        It calculates the RMSD of all the structures belonging to each
        cluster.

        Parameters
        ----------
        clusters : a numpy.array object
            The array of cluster labels that were obtained
        coordinates : numpy.array
            The array of coordinates that have been clustered

        Returns
        -------
        rmsd_per_cluster : dict[int, float]
            The mean RMSD of each cluster
        """
        from collections import defaultdict
        import numpy as np
        from pele_platform.analysis.clustering import Clustering

        coordinates = Clustering.fix_coordinates_shape(coordinates)

        if len(clusters) != len(coordinates):
            print("Warning: coordinates array has a wrong size. " +
                  "The RMSD analysis will be skipped. It will be " +
                  "skipped. Expected size: {}".format(len(clusters)))
            return

        # Split conformations by cluster
        conformations_per_cluster = defaultdict(list)
        for cluster, conformation in zip(clusters, coordinates):
            conformations_per_cluster[cluster].append(conformation)

        # Convert lists to numpy arrays
        for cluster, conformations in conformations_per_cluster.items():
            conformations_per_cluster[cluster] = np.array(conformations)

        # Calculate centroids of each cluster
        centroid_per_cluster = {}
        for cluster, conformations in conformations_per_cluster.items():
            centroid_per_cluster[cluster] = \
                np.mean(conformations, axis=0).reshape(-1, 3)

        # Calculate mean RMSD of each cluster with respect to their centroid
        rmsd_per_cluster = {}

        for cluster, conformations in conformations_per_cluster.items():
            rmsds = []
            for conformation in conformations:
                diff = \
                    conformation.reshape(-1, 3) - centroid_per_cluster[cluster]
                norm_factor = len(centroid_per_cluster[cluster])
                rmsds.append(np.sqrt((diff * diff).sum() / norm_factor))

            rmsd_per_cluster[cluster] = np.mean(rmsds)

        return rmsd_per_cluster

    def _save_cluster_representatives(self, clusters, dataframe, path,
                                      representatives_criterion):
        """
        It saves the resulting clusters to disk. The selection of the
        representative structures is based on the total energy. The
        structure chosen as the representative will be the closer
        to the 5th percentile of total energy.

        Parameters
        ----------
        clusters : a numpy.array object
            The array of cluster labels that were obtained
        dataframe : a pandas.dataframe object
            The dataframe containing the PELE reports information that
            follows the same ordering as the array of clusters
        path : str
            The path where the clusters will be saved at
        representatives_criterion : str
            Criterion to select cluster representative structures.
            One of ["total_25_percentile", "total_5_percentile",
            "total_mean", "total_min", "interaction_25_percentile",
            "interaction_5_percentile", "interaction_mean",
            "interaction_min"]
        """
        import os
        from collections import defaultdict
        import numpy as np

        from pele_platform.Utilities.Helpers import get_suffix
        from pele_platform.Utilities.Helpers.bestStructs import (
            extract_snapshot_from_pdb,
            extract_snapshot_from_xtc)
        from pele_platform.analysis.clustering import get_cluster_label
        from pele_platform.constants.constants import \
            cluster_representatives_criterion

        # Get metric to be used in the cluster representatives selection
        representatives_criterion = representatives_criterion.lower()
        if representatives_criterion not in cluster_representatives_criterion:
            raise ValueError('Invalid cluster_representatives_criterion'
                             '\'{}\''.format(representatives_criterion) +
                             '. It must be one of ' +
                             '{}'.format(
                                 cluster_representatives_criterion.keys()))

        if representatives_criterion.startswith('total'):
            metric = 'currentEnergy'
        else:
            metric = dataframe.columns[self.be_column - 1]

        # Get Binding Energy per cluster
        metrics = list(dataframe[metric])
        metrics_per_cluster = defaultdict(list)
        for cluster, value in zip(clusters, metrics):
            # Skip outliers such as clusters with label -1
            if cluster < 0:
                continue
            metrics_per_cluster[cluster].append(value)

        golden_values_per_cluster = {}
        if '_25_percentile' in representatives_criterion:
            for cluster, metrics_array in metrics_per_cluster.items():
                golden_values_per_cluster[cluster] = \
                    np.percentile(metrics_array, 25)
        elif '_5_percentile' in representatives_criterion:
            for cluster, metrics_array in metrics_per_cluster.items():
                golden_values_per_cluster[cluster] = \
                    np.percentile(metrics_array, 5)
        elif '_min' in representatives_criterion:
            for cluster, metrics_array in metrics_per_cluster.items():
                golden_values_per_cluster[cluster] = np.min(metrics_array)
        else:
            for cluster, metrics_array in metrics_per_cluster.items():
                golden_values_per_cluster[cluster] = np.mean(metrics_array)

        representative_structures = {}
        lowest_energetic_diff = {}
        trajectories = list(dataframe["trajectory"])
        steps = list(dataframe["numberOfAcceptedPeleSteps"])
        for cluster, metric, trajectory, step in zip(clusters, metrics,
                                                     trajectories, steps):
            # Skip outliers such as clusters with label -1
            if cluster < 0:
                continue

            energetic_diff = \
                np.abs(golden_values_per_cluster[cluster] - metric)
            if cluster not in representative_structures:
                representative_structures[cluster] = [trajectory, step]
                lowest_energetic_diff[cluster] = energetic_diff

            elif lowest_energetic_diff[cluster] > energetic_diff:
                representative_structures[cluster] = [trajectory, step]
                lowest_energetic_diff[cluster] = energetic_diff

        for cluster, [trajectory, step] in representative_structures.items():
            if not self.topology:
                try:
                    label = get_cluster_label(cluster)
                    extract_snapshot_from_pdb(
                        path=trajectory,
                        f_id=get_suffix(os.path.splitext(trajectory)[0]),
                        output=path,
                        topology=self.topology,
                        step=step,
                        out_freq=1,
                        f_out="cluster_{}.pdb".format(label))
                except UnicodeDecodeError:
                    raise Exception("XTC output being treated as PDB. " +
                                    "Please specify XTC with the next " +
                                    "flag. traj: 'trajectory_name.xtc' " +
                                    "in the input.yaml")
            else:
                label = get_cluster_label(cluster)
                extract_snapshot_from_xtc(
                    path=trajectory,
                    f_id=get_suffix(os.path.splitext(trajectory)[0]),
                    output=path,
                    topology=self.topology,
                    step=step,
                    out_freq=1,
                    f_out="cluster_{}.pdb".format(label))

        self._save_top_selections(representative_structures, path, dataframe)

    def _save_top_selections(self, representative_structures, path, dataframe):
        """
        It saves trajectory information about cluster representatives to
        a CSV file, then joins that data with metric from cluster summary
        dataframe (containing energies, percentiles, etc.).

        Parameters
        ----------
        representative_structures : dict[int, tuple[str, int]]
            Dictionary where cluster ID is the key and value is a list
            with [trajectory, step] of each cluster
        path : str
            The path where the CSV file will be saved at
        dataframe : pandas.DataFrame
            Dataframe with data on energies, SASA, etc. for each pose
        """
        import os
        from collections import defaultdict
        import pandas as pd
        from pele_platform.analysis.clustering import get_cluster_label

        # Gather information about each representative structure
        cluster_ids = []
        epochs = []
        trajectories = []
        steps = []
        labels = []
        for cluster_id, (trajectory, step) \
                in representative_structures.items():
            cluster_ids.append(str(cluster_id))
            trajectories.append(trajectory)
            steps.append(step)
            labels.append(get_cluster_label(cluster_id))
            epoch = os.path.basename(os.path.dirname(trajectory))
            if epoch.isdigit():
                epochs.append(epoch)
            else:
                epochs.append('-')

        # Gather metrics for each representative structure
        metrics = self._data_handler.get_metrics()
        skip = False
        metric_values = defaultdict(list)
        for metric in metrics:
            if skip:
                metric_values = {}
                break
            for cluster_id, (trajectory, step) \
                    in representative_structures.items():
                filtered_df = dataframe[dataframe['trajectory'] == trajectory]
                filtered_df = \
                    filtered_df[filtered_df['numberOfAcceptedPeleSteps'] ==
                                step]
                if len(filtered_df) != 1:
                    print('Unable to find metric \'{}\' '.format(metric) +
                          'for representative structure: ' +
                          '{}, step {}'.format(trajectory, step))
                    skip = True
                    break
                metric_values[metric].append(float(filtered_df[metric]))

        # Build dataframe
        representatives_data = pd.DataFrame({"Cluster": cluster_ids,
                                             "Cluster label": labels,
                                             "epoch": epochs,
                                             "trajectory": trajectories,
                                             "Step": steps})
        for metric, values in metric_values.items():
            representatives_data[metric] = values

        # Save csv file
        file_name = os.path.join(path, "top_selections.csv")
        representatives_data.to_csv(file_name, index=False)

    @staticmethod
    def _check_existing_directory(path):
        """
        It checks if the results folder exists and enumerates a new folder
        name to avoid overwriting the analysis.

        Parameters
        ----------
        path : str
            Path to analysis working folder

        Returns
        -------
            New working folder for analysis, if 'results' already exists,
            otherwise returns the same folder
        """
        import os

        # If current path does not exist, we are done
        if not os.path.exists(path):
            return path

        # Otherwise, suggest new path
        dir_name = os.path.dirname(path)
        folder_name = os.path.basename(path)
        chunks = folder_name.split('_')
        last_chunk = chunks[-1]

        # If last chunk is digit, enumerate starting from it
        if last_chunk.isdigit():
            new_id = int(last_chunk) + 1
            folder_name = '_'.join(chunks[:-1])
        else:
            new_id = 1
            folder_name = '_'.join(chunks)

        # Add new id to folder name
        new_folder_name = folder_name + '_' + str(new_id)

        # Concatenate old directory with new folder name
        new_path = os.path.join(dir_name, new_folder_name)

        # Iterate until finding a non existing path
        while os.path.exists(new_path):
            new_id += 1
            new_folder_name = folder_name + '_' + str(new_id)
            new_path = os.path.join(dir_name, new_folder_name)

        return new_path
