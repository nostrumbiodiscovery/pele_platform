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

        plots_folder = os.path.join(path, "plots")
        top_poses_folder = os.path.join(path, "top_poses")
        clusters_folder = os.path.join(path, "clusters")

        if not os.path.exists(plots_folder):
            os.makedirs(plots_folder)
        if not os.path.exists(top_poses_folder):
            os.makedirs(top_poses_folder)
        if not os.path.exists(clusters_folder):
            os.makedirs(clusters_folder)

        self.generate_plots(plots_folder)
        self.generate_top_poses(top_poses_folder)
        self.generate_clusters(clusters_folder, clustering_type)

    def generate_plots(self, path):
        """
        It generates the plots.

        Parameters
        ----------
        path : str
            The path where the plots will be saved
        """
        from pele_platform.analysis import Plotter

        # Get dataframe
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
        if 'Binding Energy' in metrics:
            t_energy = 'currentEnergy'
            i_energy = 'Binding Energy'
        else:
            t_energy = 'currentEnergy'
            i_energy = None

        # Iterate over all the metrics found in the reports
        for metric in metrics:
            # Skip metric if it belongs to the Total energy or the
            # Interaction energy
            if metric == 'currentEnergy' or metric == 'Binding Energy':
                continue

            if i_energy is not None:
                plotter.plot_two_metrics(t_energy, i_energy, metric,
                                         output_folder=path)
                plotter.plot_two_metrics(i_energy, metric,
                                         output_folder=path)
            else:
                plotter.plot_two_metrics(t_energy, metric,
                                         output_folder=path)

    def generate_top_poses(self, path, n_poses=100):
        """
        It selects and saves the top poses.

        Parameters
        ----------
        path : str
            The path where the plots will be saved
        n_poses : int
            The number of top poses to retrieve
        """
        # Get metrics and locate Interaction energy
        metrics = self._data_handler.get_metrics()

        if 'Binding Energy' in metrics:
            metric = 'Binding Energy'
        else:
            metric = 'currentEnergy'

        self.parameters.logger.info("Retrieve 100 Best Poses")

        top_poses = self._data_handler.get_top_entries(metric, n_poses)
        self._extract_poses(top_poses, metric, path)

    def generate_clusters(self, path, clustering_type, *args, **kargs):
        """
        It generates the structural clustering of ligand poses.

        Parameters
        ----------
        path : str
            The path where the plots will be saved
        clustering_type : str
            The clustering method that will be used to generate the
            clusters
        args : list
            The list of arguments to pass to the cluster class
        kargs : dict
            The dictionary of keyword arguments to pass to the cluster class
        """
        if not self.parameters.topology:
            coordinates = self._data_handler.extract_raw_coords(
                self.parameters.residue, filter=True,
                threshold=None, remove_hydrogen=True,
                n_proc=self.parameters.cpus)
        else:
            coordinates = self._data_handler.extract_coords(
                self.parameters.residue, self.parameters.topology,
                filter=True, threshold=None, remove_hydrogen=True)

        print(coordinates)
        print(coordinates.shape)

        self.parameters.logger.info(f"Retrieve {nclusts} best cluster poses")

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





