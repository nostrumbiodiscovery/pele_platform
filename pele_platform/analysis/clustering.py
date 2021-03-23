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
    def get_clusters(self, coordinates):
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
            now have the following shape: [M, N, 3], where M is the
            total number of models that have been sampled with PELE and
            N is the total number of atoms belonging to the residue that
            is being analyzed
        """
        try:
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

        reshaped_coordinates = coordinates.reshape(-1, n_atoms * n_dimensions)
        return reshaped_coordinates


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

    def get_clusters(self, coordinates):
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

        Returns
        -------
        clusters : numpy.array
            The array of cluster labels assigned to each conformer from
            the supplied array
        """
        from sklearn.mixture import GaussianMixture

        coordinates = Clustering.fix_coordinates_shape(coordinates)

        clustering_method = GaussianMixture(n_components=self._n_clusters,
                                            covariance_type="full")
        clusters = clustering_method.fit_predict(coordinates)

        return clusters


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

    def get_clusters(self, coordinates):
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

        Returns
        -------
        clusters : numpy.array
            The array of cluster labels assigned to each conformer from
            the supplied array
        """
        from hdbscan import HDBSCAN
        coordinates = Clustering.fix_coordinates_shape(coordinates)
        clustering_method = HDBSCAN(cluster_selection_epsilon=self._bandwidth)
        clusters = clustering_method.fit_predict(coordinates)

        return clusters


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

    def get_clusters(self, coordinates):
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

        Returns
        -------
        clusters : numpy.array
            The array of cluster labels assigned to each conformer from
            the supplied array
        """
        from sklearn.cluster import MeanShift

        coordinates = Clustering.fix_coordinates_shape(coordinates)

        clustering_method = MeanShift(bandwidth=self._bandwidth,
                                      cluster_all=True,
                                      max_iter=10000)
        clusters = clustering_method.fit_predict(coordinates)

        return clusters


def get_cluster_label(cluster_id):
    """
    It assigns a cluster label according to the cluster id that is
    supplied.

    It follows the criterion from above:
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
        The id of the cluster that will be used to generate the label.
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

    return cluster_label
