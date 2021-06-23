"""
This module contains classes and methods to handle data coming from PELE
trajectories.
"""


class DataHandler(object):
    """
    Main class to handle data coming from PELE trajectories.
    """
    _NON_METRIC_LABELS = {'#Task', 'Step', 'trajectory', 'epoch',
                          'numberOfAcceptedPeleSteps'}
    _TRAJECTORY_LABEL = 'trajectory'

    def __init__(self, sim_path, report_name, trajectory_name,
                 be_column=None, skip_initial_structures=True):
        """
        It initializes a DataHandler object.

        Parameters
        ----------
        sim_path : str
            The simulation path containing the output files coming from
            PELE
        report_name : str
            The name of PELE report files
        trajectory_name : str
            The name of PELE trajectory files
        be_column : int
            The column that belongs to the Interaction energy metric in
            PELE report files. Default is None
        skip_initial_structures : bool
            Whether to skip initial structures when extracting metrics
            or coordinates or not. Default is True
        """
        self._sim_path = sim_path
        self._report_name = report_name
        self._trajectory_name = trajectory_name
        self._be_column = be_column
        self._dataframe = None
        self.skip_initial_structures = skip_initial_structures

    @classmethod
    def from_parameters(cls, parameters):
        """
        It initializes a DataHandler object from a Parameters object.

        Parameters
        ----------
        parameters : a Parameters object
            The Parameters object containing the parameters that belong
            to the simulation

        Returns
        -------
        data_handler : a DataHandler object
            The DataHandler object obtained from the parameters that were
            supplied
        """
        import os

        sim_path = os.path.join(parameters.pele_dir,
                                parameters.output)
        report_name = parameters.report_name
        trajectory_name = parameters.traj_name
        be_column = parameters.be_column

        if parameters.test is not None:
            skip_initial_structures = not parameters.test
        else:
            skip_initial_structures = True

        data_handler = DataHandler(sim_path, report_name, trajectory_name,
                                   be_column, skip_initial_structures)
        return data_handler

    @classmethod
    def from_dataframe(cls, dataframe):
        """
        It initializes a DataHandler object from a reports dataframe.

        Parameters
        ----------
        dataframe : a pandas.DataFrame object
            The dataframe containing the information from PELE reports

        Returns
        -------
        data_handler : a DataHandler object
            The DataHandler object obtained from the dataframe that was
            supplied
        """
        import os

        arbitrary_trajectory_path = list(dataframe['trajectory'])[0]

        sim_path = os.path.dirname(os.path.dirname(arbitrary_trajectory_path))
        report_name = None
        trajectory_name = os.path.basename(arbitrary_trajectory_path)

        columns = list(dataframe.columns)
        try:
            be_column = columns.index('Binding Energy') + 1
        except ValueError:
            be_column = None

        data_handler = DataHandler(sim_path, report_name, trajectory_name,
                                   be_column)
        data_handler._dataframe = dataframe

        return data_handler

    def get_reports_dataframe(self, from_scratch=False):
        """
        It returns the data stored in PELE reports as a pandas dataframe.

        Parameters
        ----------
        from_scratch : bool
            If it is set to True, a new dataframe will be generated from
            scratch. Default is False

        Returns
        -------
        dataframe : a pandas.DataFrame object
            The dataframe containing the information from PELE reports
        """
        # Return the dataframe if it has been already created,
        # unless a brand new dataset is requested

        if self._dataframe is not None and not from_scratch:
            return self._dataframe

        # This will happen when de DataHandler has been initialized from a
        # dataframe
        if self._report_name is None:
            raise Exception('Dataframe cannot be generated from scratch '
                            + 'when report names are unknown')

        import os
        import glob
        import pandas as pd
        from pele_platform.Utilities.Helpers import get_suffix

        # Initialize primary variables
        epoch_dirs = glob.glob(os.path.join(self._sim_path, '[0-9]*'))
        report_prefix = self._report_name
        trajectory_prefix = \
            str(os.path.splitext(self._trajectory_name)[0]) + '_'
        trajectory_format = \
            str(os.path.splitext(self._trajectory_name)[-1])

        # Filter out non digit folders
        epochs = [os.path.basename(path) for path in epoch_dirs
                  if os.path.basename(path).isdigit()]

        # Sort epochs by number
        epochs = sorted(epochs, key=int)

        # Tweak to read a directory from standard PELE (not coming
        # from adaptive)
        if len(epochs) == 0:
            report_dirs = glob.glob(os.path.join(self._sim_path,
                                                 report_prefix + '_[0-9]*'))
            if len(report_dirs) > 0:
                epochs = ['']

        dataframe_lists = []
        for adaptive_epoch in epochs:
            folder = os.path.join(self._sim_path, str(adaptive_epoch))
            report_dirs = glob.glob(os.path.join(folder,
                                                 report_prefix + '_[0-9]*'))

            report_ids = [get_suffix(path) for path in report_dirs
                          if get_suffix(path).isdigit()]
            report_list = [os.path.join(folder, report_prefix + '_' + i)
                           for i in sorted(report_ids, key=int)]

            if len(report_dirs) == 0:
                print('Warning: no PELE reports with prefix ' +
                      '\'{}\'were found in: '.format(report_prefix) +
                      '{}'.format(folder))

            for i, report in enumerate(report_list, start=1):
                pandas_df = pd.read_csv(report, sep="    ", engine="python",
                                        index_col=False, header=0)
                pandas_df["epoch"] = adaptive_epoch
                pandas_df["trajectory"] = \
                    os.path.join(self._sim_path, adaptive_epoch,
                                 trajectory_prefix + str(i) +
                                 trajectory_format)
                dataframe_lists.append(pandas_df)

        if len(dataframe_lists) == 0:
            raise ValueError('No PELE trajectories were found in the ' +
                             'output that was supplied: ' +
                             '{}'.format(self._sim_path))

        self._dataframe = pd.concat(dataframe_lists, ignore_index=True)

        return self._dataframe

    def remove_outliers_from_dataframe(self, dataframe, threshold=None):
        """
        Given a dataframe, it removes the outliers by deleting the entries
        with highest values for the total and binding energy.

        Parameters
        ----------
        dataframe : a pandas.DataFrame object
            The dataframe containing the information from PELE reports
        threshold : float
            The ratio of high-energy entries that will be filtered out.
            Default is None and will be initialized with a threshold of
            0.02

        Returns
        -------
        dataframe : a pandas.DataFrame object
            The filtered dataframe containing the information from PELE
            reports
        """
        # Check threshold value
        if threshold is None:
            threshold = 0.02
        elif threshold > 1 or threshold < 0:
            raise ValueError('Invalid threshold value: '
                             'it must be higher than 0 and smaller than 1.')

        # Get the number of entries to remove
        cols = list(dataframe.columns)
        n_points_to_remove = int(len(dataframe[cols[0]]) * threshold)

        # Remove entries with higher total energies
        dataframe_filtered = dataframe.sort_values(
            cols[3], ascending=False).iloc[n_points_to_remove:]

        # Remove entries with higher interaction energies
        if self._be_column:
            dataframe_filtered = dataframe_filtered.sort_values(
                cols[self._be_column - 1],
                ascending=False).iloc[n_points_to_remove:]

        return dataframe_filtered

    def get_metrics(self):
        """
        It returns the labels that belong to the metrics in the dataframe.

        Returns
        -------
        metrics : list[str]
            The list of metrics that belong to the reports dataframe
        """
        # Get dataframe
        dataframe = self.get_reports_dataframe()

        # Get columns
        columns = list(dataframe.columns)

        # Filter out non metric columns
        metrics = [metric for metric in columns
                   if metric not in self._NON_METRIC_LABELS]

        return metrics

    def get_number_of_metrics(self):
        """
        It returns the number of metrics in the dataset.

        Returns
        -------
        n_metrics : int
            The number of metrics in the dataset
        """
        # Get metrics
        metrics = self.get_metrics()

        # Calculate number of metrics
        n_metrics = len(metrics)

        return n_metrics

    def get_top_entries(self, metric, n_entries, criterion='lowest'):
        """
        It returns the top entries according to the supplied parameters.

        Parameters
        ----------
        metric : str
            The metric to evaluate
        n_entries : int
            The number of entries to return
        criterion : str
            The criterion to evaluate the metrics. One of ['lowest',
            'largest']. If 'lowest' the best entries will be those with
            the lowest values. If 'largest' the entries whose values are
            the largest will be retrieved. Default is 'lowest'

        Returns
        -------
        dataframe :  a pandas.DataFrame object
            The dataframe containing the top entries that were filtered
        """
        # Check metric value
        metrics = self.get_metrics()
        if not str(metric).isdigit() not in metrics:
            raise ValueError('Invalid metric: metric name not found '
                             + 'in the reports dataframe')

        # Ensure that metric is pointing to a dataframe column
        if str(metric).isdigit():
            metric = self.get_column_name(metric)

        # Check criterion value
        if criterion not in ['lowest', 'largest']:
            raise ValueError('Invalid criterion: it must be one of '
                             + '[\'smallest\', \'largest\']')

        if criterion == 'lowest':
            return self._dataframe.nsmallest(n_entries, metric)
        else:
            return self._dataframe.nlargest(n_entries, metric)

    def get_column_name(self, column_index):
        """
        It returns the column name that corresponds to the index that is
        supplied. Take into account that the index starts at 1, not at 0.

        Parameters
        ----------
        column_index : int
            The index of the column whose name will be returned. It starts
            at 1, not at 0

        Returns
        -------
        column_name : str
            The name of the column that corresponds to the index that is
            supplied
        """
        dataframe = self.get_reports_dataframe()
        column_name = list(dataframe)[int(column_index) - 1]

        return column_name

    def extract_XTC_coords(self, residue_name, topology, water_ids=[],
                           remove_hydrogen=True, max_coordinates=6):
        """
        This method employs mdtraj to extract the coordinates that
        belong to the supplied residue from all the trajectories in the
        dataframe. It supports both PDB and XTC trajectories (although
        right now it is only used to deal with XTC)

        Parameters
        ----------
        residue_name : str
            A 3-char string that represent the residue that will be extracted
        topology : str
            Path to the PDB file representing the topology of the system
        water_ids : list[tuple[str, int]]
            The list of water ids whose coordinates will be extracted.
            Each water id is defined with a tuple that contains the PDB
            chain and the residue number corresponding to each water
            molecule to track. Default is []
        remove_hydrogen : bool
            Whether to remove all hydrogen atoms from the extracted
            coordinates array or not. Default is True
        max_coordinates : int
            The maximum number of coordinates to keep for each model. Default
            is 6

        Returns
        -------
        coordinates : numpy.array
            The array of coordinates that will be clustered. Its shape
            fulfills the following dimensions: [M, N, 3], where M is the
            total number of models that have been sampled with PELE and
            N is the total number of atoms belonging to the residue that
            is being analyzed
        water_coordinates : numpy.array
            The list of arrays with water coordinates. Their shape fulfills
            the following dimensions: [M, N, 3], where M is the total number
            of models that have been sampled with PELE and N is the total
            number of tracked water molecules (they have the same
            ordering as the supplied water_ids list)
        dataframe :  a pandas.DataFrame object
            The dataframe containing the information from PELE reports
            that matches with the array of coordinates that has been
            extracted
        """
        import mdtraj
        import numpy as np
        import pandas as pd

        from pele_platform.Utilities.Helpers import helpers

        try:
            indices_to_retrieve = self._coordinate_reduction(residue_name,
                                                             remove_hydrogen,
                                                             topology,
                                                             max_coordinates)
        except ValueError:
            return None, None, None

        # Load topology
        topology_file = topology
        topology = mdtraj.load(topology)

        # Select atom subset
        if remove_hydrogen:
            selection_str = \
                'resname == "{}" and symbol != H'.format(residue_name)
        else:
            selection_str = 'resname == "{}"'.format(residue_name)
        atom_indices = topology.top.select(selection_str)

        water_atom_indices = helpers.get_atom_indices(water_ids,
                                                      topology_file,
                                                      pdb_atom_name="OW")
        filtered_atom_indices = []

        # Apply the indices to retrieve
        for residue_index, atom_index in enumerate(atom_indices):
            if residue_index in indices_to_retrieve:
                filtered_atom_indices.append(atom_index)
        for water_index in water_atom_indices:
            filtered_atom_indices.append(water_index)

        # Get trajectories from reports dataframe
        dataframe = self.get_reports_dataframe()

        reordered_dataframe = pd.DataFrame()
        trajectories = list(set(dataframe[self._TRAJECTORY_LABEL]))

        coordinates = []
        water_coordinates = []
        for trajectory in trajectories:
            # Extract coordinates
            all_frames = mdtraj.load(trajectory, top=topology,
                                     atom_indices=filtered_atom_indices)

            residue_frames = []
            if len(water_atom_indices) != 0:
                water_frames = []
            else:
                water_frames = None

            for model_idxs in all_frames.xyz:
                residue_frames.append(model_idxs[0:len(indices_to_retrieve)] * 10)
                if water_frames is not None:
                    water_frames.append(model_idxs[-len(water_atom_indices):] * 10)

            # Reorder entries in the dataset to match with the coordinate
            # ordering
            trajectory_rows = dataframe.query(
                'trajectory=="{}"'.format(trajectory))
            trajectory_rows = trajectory_rows.sort_values(['Step'],
                                                          ascending=True)

            # Remove first entry
            if self.skip_initial_structures:
                residue_frames = residue_frames[1:]
                trajectory_rows = trajectory_rows.query('Step!=0')
                if water_frames is not None:
                    water_frames = water_frames[1:]

            # Save extracted data
            coordinates.extend(residue_frames)

            if water_frames is not None:
                water_coordinates.extend(water_frames)

            reordered_dataframe = reordered_dataframe.append(trajectory_rows)

        coordinates = np.array(coordinates)

        if len(water_ids) > 0:
            if len(water_coordinates) != len(coordinates):
                print('Warning: water coordinates could not be '
                      'extracted. Water sites will not be saved.')
                water_coordinates = None
            else:
                water_coordinates = np.array(water_coordinates)
        else:
            water_coordinates = None

        return coordinates, water_coordinates, reordered_dataframe

    # TODO this should be a static method!!!
    def extract_PDB_coords(self, residue_name, water_ids=[],
                           remove_hydrogen=True, max_coordinates=6,
                           n_proc=1):
        """
        This method extracts the the coordinates that belong to the
        supplied residue from all the trajectories in the dataframe.
        The trajectories must be written as PDB files.

        Parameters
        ----------
        residue_name : str
            A 3-char string that represent the residue that will be extracted
        water_ids : list[tuple[str, int]]
            The list of water ids whose coordinates will be extracted.
            Each water id is defined with a tuple that contains the PDB
            chain and the residue number corresponding to each water
            molecule to track. Default is []
        remove_hydrogen : bool
            Whether to remove all hydrogen atoms from the extracted
            coordinates array or not. Default is True
        max_coordinates : int
            The maximum number of coordinates to keep for each model. Default
            is 6
        n_proc : int
            The number of processors to employ to extract the coordinates.
            Default is 1, so the parallelization is deactivated

        Returns
        -------
        coordinates : numpy.array
            The array of coordinates that will be clustered. Its shape
            fulfills the following dimensions: [M, N, 3], where M is the
            total number of models that have been sampled with PELE and
            N is the total number of atoms belonging to the residue that
            is being analyzed
        water_coordinates : numpy.array
            The list of arrays with water coordinates. Their shape fulfills
            the following dimensions: [M, N, 3], where M is the total number
            of models that have been sampled with PELE and N is the total
            number of tracked water molecules (they have the same
            ordering as the supplied water_ids list)
        dataframe :  a pandas.DataFrame object
            The dataframe containing the information from PELE reports
            that matches with the array of coordinates that has been
            extracted
        """
        import pandas as pd

        no_multiprocessing = False
        try:
            from multiprocessing import Pool
            from functools import partial
        except ImportError:
            no_multiprocessing = True
        import numpy as np

        dataframe = self.get_reports_dataframe()
        trajectories = list(set(dataframe[self._TRAJECTORY_LABEL]))

        if len(trajectories) == 0:
            raise ValueError('No trajectories were found in the report ' +
                             'dataframe. Are all the variables correctly ' +
                             'set?')
        else:
            try:
                indices_to_retrieve = \
                    self._coordinate_reduction(residue_name, remove_hydrogen,
                                               trajectories[0],
                                               max_coordinates)
            except ValueError:
                return None, None, None

        if no_multiprocessing or n_proc == 1:
            residue_coordinates = []
            water_coordinates = []
            for trajectory in trajectories:
                res_coords, water_coords = \
                    self._get_coordinates_from_trajectory(
                        residue_name, remove_hydrogen, trajectory,
                        indices_to_retrieve=indices_to_retrieve,
                        water_ids=water_ids)
                residue_coordinates.append(res_coords)
                water_coordinates.append(water_coords)
        else:
            parallel_function = partial(self._get_coordinates_from_trajectory,
                                        residue_name, remove_hydrogen,
                                        indices_to_retrieve=indices_to_retrieve,
                                        water_ids=water_ids)

            with Pool(n_proc) as pool:
                out = pool.map(parallel_function, trajectories)

            residue_coordinates = [elem[0] for elem in out]
            water_coordinates = [elem[1] for elem in out]

        # Remove possible empty arrays
        coord_to_remove = []
        w_coord_to_remove = []
        traj_to_remove = []
        for idx, (coordinates_array, trajectory) \
                in enumerate(zip(residue_coordinates, trajectories)):
            if len(coordinates_array) == 0:
                coord_to_remove.append(coordinates_array)
                traj_to_remove.append(trajectory)

                if len(water_coordinates) == len(residue_coordinates):
                    w_coord_to_remove.append(water_coordinates[idx])

        for coord in coord_to_remove:
            residue_coordinates.remove(coord)

        for w_coord in w_coord_to_remove:
            water_coordinates.remove(w_coord)

        for traj in traj_to_remove:
            trajectories.remove(traj)

        # In case no water coordinates were extracted
        for water_coord in water_coordinates:
            if len(water_coord) == 0:
                if len(water_ids) != 0:
                    print('Warning: water coordinates could not be '
                          'extracted. Water sites will not be saved.')
                water_coordinates = None
                break

        # In case no coordinates were extracted
        if len(residue_coordinates) == 0:
            return None, None, None

        residue_coordinates = np.concatenate(residue_coordinates)

        if water_coordinates is not None:
            water_coordinates = np.concatenate(water_coordinates)

        # Reorder entries in the dataset to match with the coordinate
        # ordering
        reordered_dataframe = pd.DataFrame()

        for trajectory in trajectories:
            # Retrieve entries belonging to this trajectory, sorted by step
            # to match with coordinates
            trajectory_rows = dataframe.query(
                'trajectory=="{}"'.format(trajectory))
            trajectory_rows = trajectory_rows.sort_values(['Step'],
                                                          ascending=True)

            # Remove first entry, if applicable
            if self.skip_initial_structures:
                trajectory_rows = trajectory_rows.query('Step!=0')

            # Append the resulting entries to the new reordered dataframe
            reordered_dataframe = \
                reordered_dataframe.append(trajectory_rows)

        return residue_coordinates, water_coordinates, reordered_dataframe

    def _coordinate_reduction(self, residue_name, remove_hydrogen,
                              trajectory, max_coordinates):
        """
        It reduces the dimensions of the coordinates array to simplify the
        clustering.

        Given a [N, 3], where N is the total number of atoms that belong
        to the selected residue, it will reduce N until reaching the
        maximum number of coordinates that is supplied. The coordinates
        that will be kept are those that belong to atoms that are further
        away from each other, which are the ones that will be more
        meaningful to identify the different binding modes in the
        clustering.

        Parameters
        ----------
        residue_name : str
            A 3-char string that represent the residue that will be extracted
        remove_hydrogen : bool
            Whether to remove all hydrogen atoms from the extracted
            coordinates array or not. Default is True
        trajectory : str
            The trajectory to extract the coordinates from
        max_coordinates : int
            The maximum number of coordinates to keep for each model

        Returns
        -------
        residue_indices : list[int]
            The indices belonging to the residue atoms that are
            further away from each other which are the most meaningful
            ones to detect its different binding modes
        """
        import numpy as np

        # Array to store the indices to retrieve from the residue later on
        residue_indices = []

        # Extract the coordinates
        coordinates, water_coords = \
            self._get_coordinates_from_trajectory(residue_name,
                                                  remove_hydrogen,
                                                  trajectory,
                                                  only_first_model=True)

        try:
            coordinates = coordinates[0]  # There is only one model
        except IndexError:
            return ValueError('Residue {} '.format(residue_name) +
                              'not found in PDB file ' +
                              '{}'.format(trajectory))

        # Calculate the centroid
        centroid = np.mean(coordinates, axis=0)

        # Get the point that is further away from the centroid
        best_distance = None
        best_point = None
        for index, point in enumerate(coordinates):
            diff = point - centroid
            squared_distance = (diff * diff).sum()

            if best_distance is None or best_distance < squared_distance:
                best_distance = squared_distance
                best_point = index

        # Add best point's index to the array
        residue_indices.append(best_point)

        # Keep adding points until reaching the maximum number or the
        # total length of the residue
        for _ in range(0, min(len(coordinates), max_coordinates) - 1):
            min_distance_per_index = {}
            for index, point1 in enumerate(coordinates):
                # We exclude the current point if its index is already in
                # the list
                if index in residue_indices:
                    continue

                squared_distances = []
                for residue_index in residue_indices:
                    point2 = coordinates[residue_index]
                    diff = point1 - point2
                    squared_distance = (diff * diff).sum()
                    squared_distances.append(squared_distance)

                min_distance_per_index[index] = sorted(squared_distances)[0]

            index, _ = sorted(min_distance_per_index.items(),
                              key=lambda item: item[1],
                              reverse=True)[0]

            # Add new best point's index to the array
            residue_indices.append(index)
        return residue_indices

    def _get_coordinates_from_trajectory(self, residue_name, remove_hydrogen,
                                         trajectory, only_first_model=False,
                                         indices_to_retrieve=None,
                                         water_ids=[]):
        """
        Given the path of a trajectory, it returns the array of coordinates
        that belong to the chosen residue.

        The resulting array will have as many elements as models the
        trajectory has (excluding the first model which will be always
        skipped).

        This method is prone to be parallelized.

        .. todo ::
           * Output warnings with logger, not with print.

        Parameters
        ----------
        residue_name : str
            The name of the residue whose coordinates will be extracted
        remove_hydrogen : bool
            Whether to remove all hydrogen atoms from the extracted
            coordinates array or not. Default is True
        trajectory : str
            The trajectory to extract the coordinates from
        only_first_model : bool
            Whether to retrieve the coordinates of the first model in the
            trajectory or all of them. It is optional and its default
            value is False
        indices_to_retrieve : list[int]
            The indices of the residue to extract. Default is None and
            will extract all of them
        water_ids : list[tuple[str, int]]
            The list of water ids whose coordinates will be extracted.
            Each water id is defined with a tuple that contains the PDB
            chain and the residue number corresponding to each water
            molecule to track. Default is []

        Returns
        -------
        residue_coordinates : a numpy.Array object
            The resulting array of coordinates corresponding to the
            selected residue
        water_coordinates : a numpy.Array object
            The resulting array of coordinates corresponding to the
            supplied water ids
        """
        import numpy as np
        residue_coordinates = list()
        water_coordinates = list()

        # In case MODEL section is missing
        current_index = -1
        model_coords = []

        if only_first_model:
            skip_initial_structures = False
        else:
            skip_initial_structures = self.skip_initial_structures

        with open(trajectory) as f:
            inside_model = False
            current_model = 0

            for i, line in enumerate(f):
                if len(line) <= 6:
                    continue

                line_type = line[0:6]

                if line_type == "MODEL ":
                    if inside_model:
                        print('Warning: ENDMDL declaration for model ' +
                              '{} might be missing'.format(current_model))
                    inside_model = True
                    current_model += 1
                    current_index = -1
                    model_coords = []
                    water_coords = []

                if line_type == "ENDMDL":
                    if not inside_model:
                        print('Warning: MODEL declaration for model ' +
                              '{} might be missing'.format(current_model + 1))
                    inside_model = False

                    # Only add the current model coordinates if the array is
                    # not empty (to fulfill the dimensionality later on)
                    if len(model_coords) > 0:
                        residue_coordinates.append(np.array(model_coords))
                        model_coords = []
                    if len(water_coords) > 0:
                        water_coordinates.append(np.array(water_coords))
                        water_coords = []

                    # In case we are only interested in obtaining the
                    # coordinates of the first model, we are done
                    if only_first_model and current_model == 1:
                        break

                # First model will always be skipped, unless otherwise
                # established
                if current_model == 1 and skip_initial_structures:
                    continue

                if line_type == "ATOM  " or line_type == "HETATM":
                    current_residue_name = line[17:20]
                    current_chain = line[21]
                    current_residue_number = line[22:26]

                    # Look for residue coordinates
                    if current_residue_name == residue_name:
                        # Add one to current index (it initially equals -1)
                        current_index += 1

                        # In case we are interested in specific residue
                        # indices, retrieve only those
                        if (indices_to_retrieve is not None and
                                current_index not in indices_to_retrieve):
                            continue

                        # In case we have information about the element
                        # and we want to skip hydrogen atoms, do so
                        if remove_hydrogen and len(line) >= 78:
                            element = line[76:78]
                            element = element.strip()
                            element = element.strip(' ')
                            if element == 'H':
                                continue

                        try:
                            x = float(line[30:38])
                            y = float(line[38:46])
                            z = float(line[46:54])
                        except ValueError:
                            print('Warning: invalid PDB format found in ' +
                                  'line {}'.format(i) +
                                  'of trajectory {}. '.format(trajectory) +
                                  'Its coordinates will be skipped.')
                        point = np.array((x, y, z))
                        model_coords.append(point)

                    # Look for water coordinates
                    if (current_chain, int(current_residue_number)) in water_ids:
                        pdb_atom_name = line[12:16].strip()

                        # We only want the coordinates of oxygen atom
                        if pdb_atom_name == "OW":
                            try:
                                x = float(line[30:38])
                                y = float(line[38:46])
                                z = float(line[46:54])
                            except ValueError:
                                print('Warning: invalid PDB format found in ' +
                                      'line {}'.format(i) +
                                      'of trajectory {}. '.format(trajectory) +
                                      'Its coordinates will be skipped.')
                            point = np.array((x, y, z))
                            water_coords.append(point)
        # In case MODEL section was missing
        if not inside_model and len(model_coords) > 0:
            residue_coordinates.append(np.array(model_coords))

        if water_ids:
            if not inside_model and len(water_coords) > 0:
                water_coordinates.append(np.array(water_coords))

        residue_coordinates = np.array(residue_coordinates)
        water_coordinates = np.array(water_coordinates)

        """
        if water_ids:
            n_models, water_atoms, dimensions = water_coordinates.shape
            water_coordinates = water_coordinates.reshape(n_models * water_atoms, dimensions)
        """

        # When np.array.shape does not return a tuple of len 3 is because
        # its subarrays does not share the same dimensionality, so ligand
        # or water sizes are different.

        try:
            n_models_loaded, ligand_size, spatial_dimension = \
                residue_coordinates.shape
        except ValueError:
            if len(residue_coordinates) > 0:
                print('Warning: trajectory {} '.format(trajectory) +
                      'has an inconsistent ligand size throughout the ' +
                      'models. Its coordinates will be skipped.')

            # Return empty array
            return np.array(()), np.array(())

        if (n_models_loaded != current_model - 1 or spatial_dimension != 3) \
                and skip_initial_structures:
            print('Warning: unexpected dimensions found in the ' +
                  'coordinate array from trajectory {}. '.format(trajectory) +
                  'Its coordinates will be skipped.')

        try:
            n_models_loaded, water_atoms, spatial_dimension = \
                water_coordinates.shape
        except ValueError:
            if len(water_coordinates) > 0:
                print('Warning: trajectory {} '.format(trajectory) +
                      'has an inconsistent water size throughout the ' +
                      'models. Its coordinates will be skipped. ')

            # Return empty water coordinates array
            return residue_coordinates, np.array(())

        if (n_models_loaded != current_model - 1
            or water_atoms != len(water_ids)
                or spatial_dimension != 3) and skip_initial_structures:
            print('Warning: unexpected dimensions found in the ' +
                  'water coordinate array from trajectory ' +
                  '{}. '.format(trajectory) +
                  'Its coordinates will be skipped.')

            # Return empty array
            return residue_coordinates, np.array(())

        return residue_coordinates, water_coordinates
