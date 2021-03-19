"""
This module manages the analysis toolkit of the platform.
"""

__all__ = ["Analysis"]


class Analysis(object):
    """
    General class to manage all analysis operations.
    """

    _EPOCH_LABEL = 'epoch'
    _TRAJECTORY_LABEL = 'trajectory'
    _STEP_LABEL = 'numberOfAcceptedPeleSteps'

    def __init__(self, parameters):
        """
        It initializes an Analysis instance which it depends on
        the general Parameters class of the PELE Platform.

        Parameters
        ----------
        parameters : a Parameters object
            The Parameters object containing the parameters that belong
            to the simulation
        """
        from pele_platform.analysis import DataHandler
        self._parameters = parameters
        self._data_handler = DataHandler.from_parameters(parameters)
        self._dataframe = self._data_handler.get_reports_dataframe()

    @property
    def parameters(self):
        """
        It returns the Parameters object to analyze.

        Returns
        -------
        parameters : a Parameters object
            The Parameters object containing the parameters that belong
            to the simulation
        """
        return self._parameters

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

    def generate(self, path, clustering_type):
        """
        It runs the full analysis workflow (plots, top poses and clusters)
        and saves the results in the supplied path.

        Parameters
        ----------
        path : str
            The path where the analysis results will be saved
        clustering_type : str
            The clustering method that will be used to generate the
            clusters
        """
        import os

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
        self._dataframe.to_csv(summary_file, index=False)

        # Generate analysis results
        self.generate_plots(plots_folder)
        best_metrics = self.generate_top_poses(top_poses_folder)
        self.generate_clusters(clusters_folder, clustering_type)
        self.generate_report(plots_folder, top_poses_folder,
                             clusters_folder, best_metrics,
                             report_file)

    def generate_plots(self, path, existing_dataframe=None, colors=None):
        """
        It generates the plots.

        Parameters
        ----------
        path : str
            The path where the plots will be saved
        existing_dataframe : pandas.Dataframe
            Dataframe with data to plot.
        colors : list
            List of cluster IDs for colour mapping.
        """
        from pele_platform.analysis import Plotter

        # Get dataframe
        if existing_dataframe is None:
            dataframe = self.get_dataframe(filter=True)
        else:
            dataframe = existing_dataframe

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
        if 'Binding Energy' in metrics:
            t_energy = 'currentEnergy'
            i_energy = 'Binding Energy'
        else:
            t_energy = 'currentEnergy'
            i_energy = None

        # The minimum value for the limit column is 4, since previous
        # columns in PELE report files does not contain any metric
        if self.parameters.limit_column > 4:
            limit_column = self.parameters.limit_column - 4
        else:
            limit_column = 0

        # Iterate over all the metrics found in the reports
        for metric in metrics[limit_column:]:

            if i_energy is not None:
                plotter.plot_two_metrics(t_energy, i_energy, metric,
                                         output_folder=path)
                plotter.plot_two_metrics(metric, i_energy,
                                         output_folder=path)
            else:
                plotter.plot_two_metrics(metric, t_energy,
                                         output_folder=path)

            if self.parameters.kde:
                plotter.plot_kde(metric, i_energy, output_folder=path,
                                 kde_structs=self.parameters.kde_structs)

    def generate_top_poses(self, path, n_poses=100):
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

        if 'Binding Energy' in metrics:
            metric = 'Binding Energy'
        else:
            metric = 'currentEnergy'

        self.parameters.logger.info("Retrieve 100 Best Poses")

        top_poses = self._data_handler.get_top_entries(metric, n_poses)
        best_metrics = self._extract_poses(top_poses, metric, path)

        return best_metrics

    def generate_clusters(self, path, clustering_type):
        """
        It generates the structural clustering of ligand poses.

        Parameters
        ----------
        path : str
            The path where the clusters will be saved
        clustering_type : str
            The clustering method that will be used to generate the
            clusters
        """
        import os
        from string import ascii_uppercase
        from pele_platform.analysis import (GaussianMixtureClustering,
                                            HDBSCANClustering,
                                            MeanShiftClustering)

        self.parameters.logger.info(f"Extract coordinates for clustering")

        if clustering_type.lower() == 'gaussianmixture':
            clustering = \
                GaussianMixtureClustering(self.parameters.analysis_nclust)
            max_coordinates = 10
        elif clustering_type.lower() == 'hdbscan':
            clustering = HDBSCANClustering(self.parameters.bandwidth)
            max_coordinates = 10
        elif clustering_type.lower() == 'meanshift':
            clustering = MeanShiftClustering(self.parameters.bandwidth)
            max_coordinates = 5
        else:
            raise ValueError('Invalid clustering type: '
                             '\'{}\'. '.format(clustering_type) +
                             'It should be one of [\'GaussianMixture\', ' +
                             '\'HDBSCAN\', \'MeanShift\']')

        if not self.parameters.topology:
            coordinates, dataframe = \
                self._data_handler.extract_raw_coords(
                    self.parameters.residue, remove_hydrogen=True,
                    n_proc=self.parameters.cpus)
        else:
            coordinates, dataframe = self._data_handler.extract_coords(
                self.parameters.residue, self.parameters.topology,
                remove_hydrogen=True, max_coordinates=max_coordinates)

        if coordinates is None or dataframe is None:
            self.parameters.logger.info(f"Coordinate extraction failed, " +
                                        f"clustering analysis is skipped")
            return

        if len(coordinates) < 2:
            self.parameters.logger.info(f"Not enough coordinates, " +
                                        f"clustering analysis is skipped")
            return

        self.parameters.logger.info(f"Retrieve best cluster poses")

        if clustering_type.lower() == 'gaussianmixture':
            clustering = \
                GaussianMixtureClustering(self.parameters.analysis_nclust)
        elif clustering_type.lower() == 'hdbscan':
            clustering = HDBSCANClustering(self.parameters.bandwidth)
        elif clustering_type.lower() == 'meanshift':
            clustering = MeanShiftClustering(self.parameters.bandwidth)
        else:
            raise ValueError('Invalid clustering type: '
                             '\'{}\'. '.format(clustering_type) +
                             'It should be one of [\'GaussianMixture\', ' +
                             '\'HDBSCAN\', \'MeanShift\']')

        clusters = clustering.get_clusters(coordinates)

        rmsd_per_cluster = self._calculate_cluster_rmsds(clusters, coordinates)
        cluster_summary = self._analyze_clusters(clusters, dataframe,
                                                 rmsd_per_cluster, path)

        if len(cluster_summary) == 0:
            self.parameters.logger.info(f"No clusters could be obtained, " +
                                        f"clustering analysis is skipped")
            return

        cluster_subset, cluster_reindex_map = \
            self._select_top_clusters(clusters, cluster_summary)

        # Save cluster summary to file with information about selected labels
        cluster_summary['Selected labels'] = \
            [ascii_uppercase[cluster_reindex_map[cluster]]
             if cluster in cluster_reindex_map else '-'
             for cluster in cluster_summary['Cluster']]
        cluster_summary.to_csv(os.path.join(path, 'info.csv'), index=False)

        self._plot_clusters(cluster_subset, dataframe, cluster_summary, path)
        self._save_clusters(cluster_subset, dataframe, path)

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

        self.parameters.logger.info(
            "PDF summary report successfully written to: {}".format(report))

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
        from pele_platform.Utilities.Helpers.bestStructs \
            import extract_snapshot_from_pdb, extract_snapshot_from_xtc

        values = dataframe[metric].tolist()
        paths = dataframe[self._TRAJECTORY_LABEL].tolist()
        epochs = dataframe[self._EPOCH_LABEL].tolist()
        file_ids = [get_suffix(os.path.splitext(trajectory)[0])
                    for trajectory in paths]
        steps = list(map(int, dataframe[self._STEP_LABEL].tolist()))

        # TODO which is the purpose of this hardcoded distance_key?
        distance_key = 'distance0.5'
        if distance_key in dataframe.columns:
            dist_values = dataframe[distance_key].tolist()
            filename_template = "{}.{}.{}_BindEner{:.2f}_AtomDist{:.2f}.pdb"
            file_names = [
                filename_template.format(epoch, report, step, value, dist)
                for epoch, step, report, value, dist in zip(epochs,
                                                            steps, file_ids,
                                                            values,
                                                            dist_values)]
        else:
            filename_template = "{}.{}.{}_BindEner{:.2f}.pdb"
            file_names = [
                filename_template.format(epoch, report, step, value)
                for epoch, step, report, value in zip(epochs, steps,
                                                      file_ids, values)]

        # Read trajectory and output snapshot
        for f_id, f_out, step, path in zip(file_ids, file_names, steps, paths):
            if not self.parameters.topology:
                try:
                    extract_snapshot_from_pdb(path=path, f_id=f_id,
                                              output=output_path,
                                              topology=self.parameters.topology,
                                              step=step, out_freq=1,
                                              f_out=f_out,
                                              logger=self.parameters.logger)
                except UnicodeDecodeError:
                    raise Exception('XTC output being treated as PDB. '
                                    + 'Please specify XTC with the next '
                                    + 'flag. traj: \'trajectory_name.xtc\' '
                                    + 'in the input.yaml')
            else:
                extract_snapshot_from_xtc(path=path, f_id=f_id,
                                          output=output_path,
                                          topology=self.parameters.topology,
                                          step=step, out_freq=1,
                                          f_out=f_out,
                                          logger=self.parameters.logger)

        return values

    def _analyze_clusters(self, clusters, dataframe, rmsd_per_cluster,
                          path):
        """
        It analyzes the clusters and generates a summary with all
        the calculated descriptors. It also generates some plots
        with information about the clusters.

        .. todo ::
           * Generate box plots
           * Generate scatter plots colored by cluster
           * Calculate mean RMSD for each cluster

        Parameters
        ----------
        clusters : a numpy.array object
            The array of cluster labels that were obtained
        dataframe : a pandas.dataframe object
            The dataframe containing the PELE reports information that
            follows the same ordering as the array of clusters
        rmsd_per_cluster : dict[int, float]
            The mean RMSD of each cluster
        path : str
            The path where the output files will be saved at

        Returns
        -------
        cluster_summary : a pandas.dataframe object
            The dataframe containing summary of all clusters that were
            analyzed
        """
        import os
        from collections import defaultdict
        import pandas as pd
        import numpy as np
        from matplotlib import pyplot as plt

        metrics = self._data_handler.get_metrics()

        clusters_population = defaultdict(int)
        for cluster in clusters:
            clusters_population[cluster] += 1

        summary = pd.DataFrame([(cluster, population / len(clusters),
                                 rmsd_per_cluster[cluster])
                                for cluster, population
                                in clusters_population.items()
                                if not cluster < 0],
                               columns=['Cluster', 'Population', 'MeanRMSD'])

        # Plot Mean RMSD per cluster
        fig, ax = plt.subplots()
        ax.scatter([cluster
                    for cluster in sorted(clusters_population.keys())
                    if cluster != -1],
                   [rmsd_per_cluster[cluster]
                    for cluster in sorted(clusters_population.keys())
                    if cluster != -1],
                   s=[clusters_population[cluster] / len(clusters) * 200
                      for cluster in sorted(clusters_population.keys())
                      if cluster != -1],
                   label='Cluster population')
        ax.set_xlabel('Cluster label')
        ax.set_ylabel('Mean RMSD (Å)')
        plt.legend()
        plot_filename = os.path.join(path, "clusters_meanRMSD.png")
        plt.savefig(plot_filename)

        # Generate descriptors and boxplots for all reported metrics
        descriptors = defaultdict(dict)
        for metric in metrics:
            values_per_cluster = defaultdict(list)
            values = list(dataframe[metric])

            if len(clusters) != len(values):
                logger = self.parameters.logger
                logger.info('Warning: metric \'{}\' '.format(metric) +
                            'array has a wrong size. It will be skipped ' +
                            'from the clustering analysis. Expected size: ' +
                            '{}'.format(len(clusters)))
                continue

            # Arrange metrics per cluster
            for cluster, value in zip(clusters, values):
                # Skip outliers
                if cluster < 0:
                    continue
                values_per_cluster[cluster].append(value)

            # Calculate descriptors
            for cluster in values_per_cluster:
                descriptors['{} min'.format(metric)][cluster] = \
                    np.min(values_per_cluster[cluster])
                descriptors['{} 5-percentile'.format(metric)][cluster] = \
                    np.percentile(values_per_cluster[cluster], 5)
                descriptors['{} 25-percentile'.format(metric)][cluster] = \
                    np.percentile(values_per_cluster[cluster], 25)
                descriptors['{} mean'.format(metric)][cluster] = \
                    np.mean(values_per_cluster[cluster])
                descriptors['{} 75-percentile'.format(metric)][cluster] = \
                    np.percentile(values_per_cluster[cluster], 75)
                descriptors['{} 95-percentile'.format(metric)][cluster] = \
                    np.percentile(values_per_cluster[cluster], 95)
                descriptors['{} max'.format(metric)][cluster] = \
                    np.max(values_per_cluster[cluster])

            # Generate boxplots
            try:
                fig, ax = plt.subplots()

                ax.boxplot([values_per_cluster[cluster]
                            for cluster in sorted(values_per_cluster.keys())])

                ax.set_ylabel(metric)
                ax.set_xlabel("Cluster number")

                boxplot_filename = \
                    os.path.join(path,
                                 "clusters_{}_boxplot.png".format(metric))
                plt.savefig(boxplot_filename)
            except IndexError:
                logger = self.parameters.logger
                logger.info("Samples too disperse to produce a cluster " +
                            "for metric {}".format(metric))

        # Add descriptors to summary dataframe
        for label, values_per_cluster in descriptors.items():
            summary[label] = [values_per_cluster[cluster]
                              for cluster in summary['Cluster']
                              if not cluster < 0]

        return summary

    def _select_top_clusters(self, clusters, cluster_summary,
                             max_clusters_to_select=8,
                             min_population_to_select=0.01):
        """
        It selects the top clusters based on their binding energy. If
        this metric is not available, the cluster population will be
        used instead.

        .. todo ::
           * The user should be able to choose max_clusters_to_select
           * The user should be able to choose min_population_to_select
           * The user might also be able to choose the metric to use
             in order to select the top clusters

        Parameters
        ----------
        clusters : a numpy.array object
            The array of cluster labels that were obtained
        cluster_summary : a pandas.dataframe object
            The dataframe containing summary of all clusters that were
            analyzed
        max_clusters_to_select : int
            The maximum number of clusters to select as top. Default is 8
        min_population_to_select : float
            The minimum population the clusters must have in order to
            be selected. Default is 0.01, i.e. a population of 1%

        Returns
        -------
        cluster_subset : a numpy.array object
            The array of cluster after the selection. Those clusters
            that were not selected are now labeled with a -1
        cluster_reindex_map : dict[int, int]
            It pairs the old cluster label (dictionary key) with the new
            one after the filtering (dictionary value)
        """
        metrics = list(cluster_summary.columns)

        if 'Binding Energy 25-percentile' in metrics:
            metric = 'Binding Energy 25-percentile'
        else:
            metric = 'Population'

        filtered_cluster_summary = \
            cluster_summary[cluster_summary['Population'] >=
                            min_population_to_select]
        filtered_cluster_summary = \
            filtered_cluster_summary.nsmallest(max_clusters_to_select,
                                               metric)
        top_clusters = list(filtered_cluster_summary['Cluster'])

        cluster_reindex_map = {}
        for index, cluster in enumerate(sorted(top_clusters)):
            cluster_reindex_map[cluster] = index

        cluster_subset = []
        for cluster in clusters:
            if cluster in top_clusters:
                cluster_subset.append(cluster_reindex_map[cluster])
            else:
                cluster_subset.append(-1)

        return cluster_subset, cluster_reindex_map

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

        if 'Binding Energy' in metrics:
            energy = 'Binding Energy'
        else:
            energy = 'currentEnergy'

        # The minimum value for the limit column is 4, since previous
        # columns in PELE report files does not contain any metric
        if self.parameters.limit_column > 4:
            limit_column = self.parameters.limit_column - 4
        else:
            limit_column = 0

        # Iterate over all the metrics found in the reports
        for metric in metrics[limit_column:]:
            plotter.plot_clusters(metric, energy,
                                  output_folder=path,
                                  clusters=clusters)

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
            logger = self.parameters.logger
            logger.info('Warning: coordinates array has a wrong size. ' +
                        'The RMSD analysis will be skipped. It will be ' +
                        'skipped. Expected size: {}'.format(len(clusters)))
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
                diff = conformation.reshape(-1, 3) - centroid_per_cluster[cluster]
                norm_factor = len(centroid_per_cluster[cluster])
                rmsds.append(np.sqrt((diff * diff).sum() / norm_factor))
            rmsd_per_cluster[cluster] = np.mean(rmsds)

        return rmsd_per_cluster

    def _save_clusters(self, clusters, dataframe, path):
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
        """
        import os
        from collections import defaultdict
        import numpy as np
        from string import ascii_uppercase
        from pele_platform.Utilities.Helpers import get_suffix
        from pele_platform.Utilities.Helpers.bestStructs \
            import extract_snapshot_from_pdb, extract_snapshot_from_xtc

        energies = list(dataframe['currentEnergy'])
        energies_per_cluster = defaultdict(list)
        for cluster, energy in zip(clusters, energies):
            # Skip outliers such as clusters with label -1
            if cluster < 0:
                continue
            energies_per_cluster[cluster].append(energy)

        percentiles_per_cluster = {}
        for cluster, energies_array in energies_per_cluster.items():
            percentiles_per_cluster[cluster] = np.percentile(energies_array, 5)

        representative_structures = {}
        lowest_energetic_diff = {}
        trajectories = list(dataframe['trajectory'])
        steps = list(dataframe['numberOfAcceptedPeleSteps'])
        for cluster, energy, trajectory, step in zip(clusters, energies,
                                                     trajectories, steps):
            # Skip outliers such as clusters with label -1
            if cluster < 0:
                continue

            energetic_diff = np.abs(percentiles_per_cluster[cluster] - energy)
            if cluster not in representative_structures:
                representative_structures[cluster] = (trajectory, step)
                lowest_energetic_diff[cluster] = energetic_diff

            elif lowest_energetic_diff[cluster] > energetic_diff:
                representative_structures[cluster] = (trajectory, step)
                lowest_energetic_diff[cluster] = energetic_diff

        for cluster, (trajectory, step) in representative_structures.items():
            if not self.parameters.topology:
                try:
                    extract_snapshot_from_pdb(
                        path=trajectory,
                        f_id=get_suffix(os.path.splitext(trajectory)[0]),
                        output=path, topology=self.parameters.topology,
                        step=step, out_freq=1,
                        f_out='cluster_{}.pdb'.format(ascii_uppercase[cluster]),
                        logger=self.parameters.logger)
                except UnicodeDecodeError:
                    raise Exception('XTC output being treated as PDB. '
                                    + 'Please specify XTC with the next '
                                    + 'flag. traj: \'trajectory_name.xtc\' '
                                    + 'in the input.yaml')
            else:
                extract_snapshot_from_xtc(
                    path=trajectory,
                    f_id=get_suffix(os.path.splitext(trajectory)[0]),
                    output=path, topology=self.parameters.topology,
                    step=step, out_freq=1,
                    f_out='cluster_{}.pdb'.format(ascii_uppercase[cluster]),
                    logger=self.parameters.logger)
