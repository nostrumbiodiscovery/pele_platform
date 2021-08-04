"""
This module contains classes and methods to cluster PELE trajectories.
"""

from abc import ABC, abstractmethod


class Clustering(ABC):
    """
    Abstract class for a Clustering method.
    """

    def __init__(self):
        """
        It initializes a Clustering object.
        """
        super().__init__()

    @abstractmethod
    def get_clusters(self, coordinates, original_df=None,
                     coordinates_df=None, csv_path=None):
        """
        It employs a clustering method to gather the supplied coordinates
        into clusters.

        Parameters
        ----------
        coordinates : numpy.array
            The array of coordinates that will be clustered. Its shape must
            fulfill the following dimensions: [M, N, 3], where M is the
            total number of models that have been sampled with PELE and
            N is the total number of atoms belonging to the residue that
            is being analyzed
        original_df : pandas.DataFrame
            Original dataframe from Analysis to be overwritten. Optional
            parameter, if not supplied the dataframe will not be updated
            with the information of clusters
        coordinates_df : pandas.DataFrame
            The filtered dataframe which was used to extract coordinates for
            clustering. Optional parameter, if not supplied the CSV
            containing the information of clusters will not be saved
        csv_path : str
            Directory where the CSV will be saved. Optional parameter,
            if not supplied the CSV containing the information of clusters
            will not be saved

        Returns
        -------
        clusters : numpy.array
            The array of cluster labels assigned to each conformer from
            the supplied array
        """
        pass

    @staticmethod
    def fix_coordinates_shape(coordinates):
        """
        Given an array of coordinates with the following dimensions:
        [M, N, 3], it reshapes it to [M, N * 3]. This is the shape
        that clustering methods require.

        Parameters
        ----------
        coordinates : numpy.array
            The array of coordinates that will be clustered. Its shape must
            fulfill the following dimensions: [M, N, 3], where M is the
            total number of models that have been sampled with PELE and
            N is the total number of atoms belonging to the residue that
            is being analyzed

        Returns
        -------
        reshaped_coordinates : numpy.array
            The reshaped array of coordinates that will be clustered. They
            now have the following shape: [M, N * 3], where M is the
            total number of models that have been sampled with PELE and
            N is the total number of atoms belonging to the residue that
            is being analyzed
        """
        try:
            if len(coordinates.shape) == 2:
                n_atoms, n_dimensions = coordinates.shape
            else:
                n_models, n_atoms, n_dimensions = coordinates.shape
            if n_dimensions != 3:
                raise ValueError
        except ValueError:
            raise ValueError('Array of coordinates have invalid ' +
                             'dimensions: {}. '.format(coordinates.shape) +
                             'Its shape must fulfill the following ' +
                             'dimensions: [M, N, 3], where M is the ' +
                             'total number of models that have been ' +
                             'sampled with PELE and N is the total ' +
                             'number of atoms belonging to the residue ' +
                             'that is being analyzed')

        if len(coordinates.shape) == 2:
            return coordinates

        reshaped_coordinates = coordinates.reshape(-1, n_atoms * n_dimensions)
        return reshaped_coordinates

    @staticmethod
    def _save_cluster_info(original_dataframe, clustering_dataframe, clusters,
                           csv_path):
        """
        Joins the original dataframe (self._dataframe from Analysis)
        with the one passed to the clustering and obtained cluster labels.
        Finally, it saves the resulting dataframe as a csv file.

        Parameters
        ----------
        original_dataframe : pandas.DataFrame
            Original dataframe from Analysis
        clustering_dataframe :
            Filtered dataframe used for coordinates extraction for the
            clustering algorithm
        clusters : np.array
            Clustering labels
        csv_path : str
            Directory where the CSV will be saved
        """
        import pandas as pd
        import os

        path = os.path.join(csv_path, "data.csv")

        clustering_dataframe.insert(len(clustering_dataframe.columns),
                                    "Cluster",
                                    [str(element) for element
                                     in clusters.tolist()])

        keys = [column for column in clustering_dataframe
                if column in original_dataframe]
        final_df = pd.merge(original_dataframe, clustering_dataframe,
                            how="left", on=keys)
        final_df = final_df.drop("#Task", axis=1).fillna("-")

        # Save csv file
        final_df.to_csv(path, index=False)


class GaussianMixtureClustering(Clustering):
    """
    Class that defines the Gaussian Mixture clustering method.
    """

    def __init__(self, n_clusters):
        """
        It initializes a GaussianMixtureClustering object.

        Parameters
        ----------
        n_clusters : int
            The number of clusters to generate
        """
        super().__init__()
        self._n_clusters = n_clusters

    def get_clusters(self, coordinates, original_df=None,
                     coordinates_df=None, csv_path=None):
        """
        It employs the Gaussian Mixture method to gather the supplied
        coordinates into clusters.

        Parameters
        ----------
        coordinates : numpy.array
            The array of coordinates that will be clustered. Its shape must
            fulfill the following dimensions: [M, N, 3], where M is the
            total number of models that have been sampled with PELE and
            N is the total number of atoms belonging to the residue that
            is being analyzed
        original_df : pandas.DataFrame
            Original dataframe from Analysis to be overwritten. Optional
            parameter, if not supplied the dataframe will not be updated
            with the information of clusters
        coordinates_df : pandas.DataFrame
            The filtered dataframe which was used to extract coordinates for
            clustering. Optional parameter, if not supplied the CSV
            containing the information of clusters will not be saved
        csv_path : str
            Directory where the CSV will be saved. Optional parameter,
            if not supplied the CSV containing the information of clusters
            will not be saved

        Returns
        -------
        clusters : numpy.array
            The array of cluster labels assigned to each conformer from
            the supplied array
        """
        from sklearn.mixture import GaussianMixture

        # Run clustering
        coordinates = Clustering.fix_coordinates_shape(coordinates)
        clustering_method = GaussianMixture(n_components=self._n_clusters, covariance_type="full")
        clusters = clustering_method.fit_predict(coordinates)

        # Save cluster information (optional)
        if (original_df is not None and coordinates_df is not None and
                csv_path is not None):
            if len(coordinates_df) > 0:
                self._save_cluster_info(original_df, coordinates_df,
                                        clusters, csv_path)

        return clusters, clustering_method


class HDBSCANClustering(Clustering):
    """
    Class that defines the HDBSCAN clustering method.
    """

    def __init__(self, bandwidth):
        """
        It initializes a HDBSCANClustering object.

        Parameters
        ----------
        bandwidth : float
            The bandwidth to employ when building clusters, it manages
            their size
        """
        super().__init__()
        self._bandwidth = bandwidth

    def get_clusters(self, coordinates, original_df=None,
                     coordinates_df=None, csv_path=None):
        """
        It employs the HDBSCAN method to gather the supplied coordinates
        into clusters.

        Parameters
        ----------
        coordinates : numpy.array
            The array of coordinates that will be clustered. Its shape must
            fulfill the following dimensions: [M, N, 3], where M is the
            total number of models that have been sampled with PELE and
            N is the total number of atoms belonging to the residue that
            is being analyzed
        original_df : pandas.DataFrame
            Original dataframe from Analysis to be overwritten. Optional
            parameter, if not supplied the dataframe will not be updated
            with the information of clusters
        coordinates_df : pandas.DataFrame
            The filtered dataframe which was used to extract coordinates for
            clustering. Optional parameter, if not supplied the CSV
            containing the information of clusters will not be saved
        csv_path : str
            Directory where the CSV will be saved. Optional parameter,
            if not supplied the CSV containing the information of clusters
            will not be saved

        Returns
        -------
        clusters : numpy.array
            The array of cluster labels assigned to each conformer from
            the supplied array
        """
        from hdbscan import HDBSCAN

        # Run clustering
        coordinates = Clustering.fix_coordinates_shape(coordinates)
        clustering_method = HDBSCAN(cluster_selection_epsilon=self._bandwidth)
        clusters = clustering_method.fit_predict(coordinates)

        # Save cluster information (optional)
        if (original_df is not None and coordinates_df is not None and
                csv_path is not None):
            if len(coordinates_df) > 0:
                self._save_cluster_info(original_df, coordinates_df,
                                        clusters, csv_path)

        return clusters, clustering_method


class MeanShiftClustering(Clustering):
    """
    Class that defines the Mean Shift clustering method.
    """

    def __init__(self, bandwidth):
        """
        It initializes a MeanShiftClustering object.

        Parameters
        ----------
        bandwidth : float
            The bandwidth to employ when building clusters, it manages
            their size
        """
        super().__init__()
        self._bandwidth = bandwidth

    def get_clusters(self, coordinates, original_df=None,
                     coordinates_df=None, csv_path=None):
        """
        It employs the Mean Shift method to gather the supplied coordinates
        into clusters.

        Parameters
        ----------
        coordinates : numpy.array
            The array of coordinates that will be clustered. Its shape must
            fulfill the following dimensions: [M, N, 3], where M is the
            total number of models that have been sampled with PELE and
            N is the total number of atoms belonging to the residue that
            is being analyzed
        original_df : pandas.DataFrame
            Original dataframe from Analysis to be overwritten. Optional
            parameter, if not supplied the dataframe will not be updated
            with the information of clusters
        coordinates_df : pandas.DataFrame
            The filtered dataframe which was used to extract coordinates for
            clustering. Optional parameter, if not supplied the CSV
            containing the information of clusters will not be saved
        csv_path : str
            Directory where the CSV will be saved. Optional parameter,
            if not supplied the CSV containing the information of clusters
            will not be saved

        Returns
        -------
        clusters : numpy.array
            The array of cluster labels assigned to each conformer from
            the supplied array
        """
        from sklearn.cluster import MeanShift

        # Run clustering
        coordinates = Clustering.fix_coordinates_shape(coordinates)
        clustering_method = MeanShift(bandwidth=self._bandwidth,
                                      cluster_all=True,
                                      max_iter=10000)
        clusters = clustering_method.fit_predict(coordinates)

        # Save cluster information (optional)
        if (original_df is not None and coordinates_df is not None and
                csv_path is not None):
            if len(coordinates_df) > 0:
                self._save_cluster_info(original_df, coordinates_df,
                                        clusters, csv_path)

        return clusters, clustering_method


def get_cluster_label(cluster_id, uppercase=True):
    """
    It assigns a cluster label according to the cluster id that is
    supplied.

    It follows the criterion from below:

    Cluster id   |   Cluster label
    0           -->  A
    1           -->  B
    2           -->  C
    25          -->  Z
    26          -->  AA
    27          -->  AB
    28          -->  AC

    Parameters
    ----------
    cluster_id : int
        The id of the cluster that will be used to generate the label
    uppercase : bool
        Whether to retrieve an uppercase label or not. Default is True

    Returns
    -------
    cluster_label : str
        The cluster label according to the supplied id and the criterion
        mentioned above
    """
    from string import ascii_uppercase

    cluster_label = ''
    current_index = cluster_id
    while current_index >= 0:
        if current_index < len(ascii_uppercase):
            cluster_label += ascii_uppercase[current_index]
        else:
            for letter in reversed(cluster_label):
                if letter != 'Z':
                    idx = ascii_uppercase.index(cluster_label[-1])
                    cluster_label = \
                        cluster_label[:-1] + ascii_uppercase[idx + 1]
                    break
            else:
                cluster_label = 'A' + cluster_label

        current_index -= 26

    if not uppercase:
        cluster_label = cluster_label.lower()

    return cluster_label
