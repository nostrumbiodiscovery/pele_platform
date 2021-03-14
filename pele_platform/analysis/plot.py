import os
import shutil
import glob
import subprocess
import mdtraj
import numpy as np
import pandas as pd
from sklearn import mixture, cluster
from multiprocessing import Pool
import AdaptivePELE.analysis.selectOnPlot as sp
import pele_platform.Utilities.Helpers.bestStructs as bs
from pele_platform.constants import constants as cs
import pele_platform.analysis.pdf_report as pr
from pele_platform.Utilities.Helpers.helpers import backup_logger

EPOCH = "epoch"
STEPS = 3
TRAJECTORY = "trajectory"


def _extract_coords(info):
    p, v, resname, topology = info
    # Most time consuming step 0.1
    traj = mdtraj.load_frame(p, v, top=topology)
    atoms_info = traj.topology.to_dataframe()[0]
    condition = atoms_info["resName"] == resname
    atom_numbers_ligand = atoms_info[condition].index.tolist()
    coords = []
    for atom_num in atom_numbers_ligand:
        try:
            coords.extend(traj.xyz[0][atom_num].tolist())
        except IndexError:
            continue
    return np.array(coords).ravel() * 10




class Plotter(object):
    """
    It handles the plots.
    """

    def __init__(self, dataframe, logger=None):
        """
        It initializes the plotter with the dataframe that contains the
        data to plot.

        Parameters
        ----------
        dataframe : a pandas.DataFrame object
            The dataframe containing the information from PELE reports
        logger : a Logger object
            The logger manager to stream out the log messages. Default is
            None
        """
        self._dataframe = dataframe
        self._logger = logger

    def plot_two_metrics(self, metric_to_x, metric_to_y, metric_to_z=None,
                         output_name=None, output_folder="."):
        """
        Given 2 or 3 metrics, it generates the scatter plot. In case that
        a 3rd metric is supplied, it will be represented as the color bar.

        Parameters
        ----------
        metric_to_x : str or int
            The metric id to plot in the X axis. It can be either a string
            with the name of the metric or an integer with the column index
        metric_to_y : str or int
            The metric id to plot in the Y axis. It can be either a string
            with the name of the metric or an integer with the column index
        metric_to_z : str or int
            The metric id to plot in the color bar. It can be either a string
            with the name of the metric or an integer with the column index.
            Default is None
        output_name : str
            The name that will be given to the resulting plot. Default is None
        output_folder : str
            The path where the plot will be saved. Default is '.', so it
            will be stored in the local directory
        """
        # Ensure that output_folder exists
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

        # Ensure that metrics are strings pointing to dataframe columns
        if str(metric_to_x).isdigit():
            metric_to_x = self._get_column_name(self._dataframe, metric_to_x)
        if str(metric_to_y).isdigit():
            metric_to_y = self._get_column_name(self._dataframe, metric_to_y)
        if metric_to_z is not None and str(metric_to_z).isdigit():
            metric_to_z = self._get_column_name(self._dataframe, metric_to_z)

        # Prepare plot name
        if output_name is None:
            if metric_to_z is not None:
                output_name = "{}_{}_{}_plot.png".format(metric_to_x,
                                                         metric_to_y,
                                                         metric_to_z)
            else:
                output_name = "{}_{}_plot.png".format(metric_to_x, metric_to_y)

        # Replace whitespaces in the output name
        output_name = os.path.join(output_folder, output_name).replace(" ", "_")

        # Generate plot with matplotlib
        from matplotlib import pyplot as plt

        fig, ax = plt.subplots()
        if metric_to_z is not None:
            scatter = ax.scatter(self._dataframe[metric_to_x],
                                 self._dataframe[metric_to_y],
                                 c=self._dataframe[metric_to_z],
                                 s=20)
            cbar = plt.colorbar(scatter)
            cbar.ax.set_ylabel(metric_to_z)
            ax.set_xlabel(metric_to_x)
            ax.set_ylabel(metric_to_y)
            plt.savefig(output_name)
            backup_logger(self._logger,
                          "Plotted {} vs {} vs {}".format(metric_to_x,
                                                          metric_to_y,
                                                          metric_to_z))
        else:
            ax.scatter(self._dataframe[metric_to_x],
                       self._dataframe[metric_to_y],
                       s=20)
            ax.set_xlabel(metric_to_x)
            ax.set_ylabel(metric_to_y)
            plt.savefig(output_name)
            backup_logger(self._logger, "Plotted {} vs {}".format(metric_to_x,
                                                                  metric_to_y))
        return output_name

    def cluster_poses(self, n_structs, metric, output, nclusts=10):
        assert self.residue, "Set residue ligand name to clusterize"
        metric = (
            metric
            if not str(metric).isdigit()
            else self._get_column_name(self.data, metric)
        )
        best_poses = self.data.nsmallest(n_structs, metric)
        clusters = self._cluster(best_poses, metric, output, nclusts)
        return clusters

    def _cluster(self, poses, metric, output, nclusts):
        # Extract metric values
        values = poses[metric].tolist()
        epochs = poses[EPOCH].tolist()
        paths = poses[TRAJECTORY].tolist()
        file_ids = [traj[:-4].split("_")[-1] for traj in paths]
        steps = self._get_column_name(self.data, STEPS)
        snapshots = poses[steps].tolist()

        # Extract coords
        input_pool = [
            [p, v, self.residue, self.topology] for p, v in zip(paths, snapshots)
        ]

        with Pool(processes=self.cpus) as pool:
            self.all_coords = pool.map(_extract_coords, input_pool)

        # Cluster
        assert self.all_coords[0][
            0
        ], "Ligand not found check the option --resname. i.e python interactive.py 5 6 7 --resname LIG"

        try:
            if self.clustering_method.lower() == "meanshift":
                self.meanshift_clustering()
            elif self.clustering_method.lower() == "dbscan":
                self.hdbscan_clustering()
            else:
                self.gmm_clustering(nclusts)
        except ValueError:
            self.labels = [1]

        n_clusters = len(set(self.labels))

        files_out = [
            "cluster{}_{}.{}.{}_BindEner{:.2f}.pdb".format(
                cluster, epoch, report, int(step), value
            )
            for epoch, step, report, value, cluster in zip(
                epochs, snapshots, file_ids, values, self.labels
            )
        ]
        all_metrics = []
        output_clusters = []
        cluster_range = range(min(self.labels), min(self.labels) + n_clusters)
        cluster_indexes = []

        for n_cluster in cluster_range:
            cluster_indexes.append(n_cluster)
            metrics = {
                value: idx
                for idx, (value, cluster) in enumerate(zip(values, self.labels))
                if n_cluster == cluster
            }
            out_freq = 1
            cluster_metrics = list(metrics.keys())
            max_idx = metrics[np.min(cluster_metrics)]
            max_traj = paths[max_idx]
            max_snapshot = snapshots[max_idx]
            output_traj = files_out[max_idx]
            input_traj = file_ids[max_idx]

            if not self.topology:
                bs.extract_snapshot_from_pdb(
                    max_traj,
                    input_traj,
                    output,
                    self.topology,
                    max_snapshot,
                    out_freq,
                    output_traj,
                    logger=self.logger,
                )
            else:
                bs.extract_snapshot_from_xtc(
                    max_traj,
                    input_traj,
                    output,
                    self.topology,
                    max_snapshot,
                    out_freq,
                    output_traj,
                    logger=self.logger,
                )
            all_metrics.append(cluster_metrics)
            output_clusters.append(os.path.join(output, output_traj))
        fig, ax = plt.subplots()
        try:
            ax.boxplot(all_metrics)
        except IndexError:
            self.logger.info("Samples to disperse to produce a cluster")
            return
        ax.set_ylabel(metric)
        ax.set_xlabel("Cluster number")
        plt.savefig(os.path.join(output, "clusters_{}_boxplot.png".format(metric.replace(" ", "_"))))

        self.write_report(output, cluster_indexes, output_clusters, all_metrics)

        return output_clusters

    @staticmethod
    def write_report(output, cluster_indexes, output_clusters, all_metrics):
        """
        Writes a report with cluster metrics, such as population, mean binding energy, etc.
        Parameters
        ----------
        output : str
            directory to save the report
        cluster_indexes : List[int]
            list of cluster indexes
        output_clusters : List[str]
            list of cluster representatives
        all_metrics : List[float]
            list of lists containing binding energies for structures belonging to each cluster
        Returns
        -------
            CSV report.
        """
        report_name = os.path.join(output, "clustering_report.csv")

        data = {"Cluster_ID": cluster_indexes,
                "Cluster_representative": output_clusters,
                "Population": [len(element) for element in all_metrics],
                "Mean_binding_energy": [np.mean(element) for element in all_metrics],
                "Max_binding_energy": [max(element) for element in all_metrics],
                "Min_binding_energy": [min(element) for element in all_metrics]}

        clusters_report = pd.DataFrame.from_dict(data)
        clusters_report.to_csv(report_name, index=False)
        return report_name

    def _get_column_name(self, df, column_digit):
        return list(df)[int(column_digit) - 1]

    def _moved_folder(self, input_file):
        if not os.path.exists(input_file):
            return True
        else:
            with open(input_file, "r") as f:
                if self.simulation_path in "".join(f.readlines()):
                    return False
                else:
                    return True

    def gmm_clustering(self, nclusts):
        """
        Performs GaussianMixture clustering on ligand coordinates, default number of clusters is 10.
        Parameters
        ----------
        nclusts : int
            Number of clusters
        Returns
        -------
            List of cluster IDs.
        """
        clf = mixture.GaussianMixture(n_components=nclusts, covariance_type="full")
        self.labels = clf.fit_predict(self.all_coords).tolist()

    def hdbscan_clustering(self):
        """
        Performs DBSCAN clustering on all ligand coordinates using user-defined epsilon (we reuse bandwidth flag).
        Returns
        -------
            List of cluster IDs.
        """
        import hdbscan
        clf = hdbscan.HDBSCAN(cluster_selection_epsilon=self.bandwidth)
        clf.fit(self.all_coords)
        self.labels = clf.labels_

    def meanshift_clustering(self):
        """
        Performs meanshift clustering on all ligand coordinates using user-defined bandwidth.
        Returns
        -------
            List of cluster IDs.
        """
        clf = cluster.MeanShift(bandwidth=self.bandwidth, cluster_all=False)
        self.labels = clf.fit_predict(self.all_coords)


def analyse_simulation(
        report_name,
        traj_name,
        simulation_path,
        residue,
        output_folder=".",
        cpus=5,
        clustering=True,
        mae=False,
        nclusts=10,
        overwrite=False,
        topology=False,
        be_column=4,
        limit_column=6,
        te_column=3,
        clustering_method="GaussianMixture",
        bandwidth=None,
        logger=None,
):
    results_folder = os.path.join(output_folder, "results")
    if os.path.exists(results_folder):
        if not overwrite:
            raise ValueError(
                "Analysis folder {} already exists. Use the option overwrite_analysis: true".format(
                    results_folder
                )
            )
        else:
            shutil.rmtree(os.path.join(output_folder, "results"))
    analysis = PostProcessor(
        report_name,
        traj_name,
        simulation_path,
        cpus,
        residue=residue,
        topology=topology,
        be_column=be_column,
        limit_column=limit_column,
        te_column=te_column,
        clustering_method=clustering_method,
        bandwidth=bandwidth,
        logger=logger,
    )

    metrics = len(list(analysis.data)) - 1  # Discard epoch as metric
    be = analysis.be_column
    total_energy = analysis.te_column
    current_metric = analysis.limit_column
    plots_folder = os.path.join(output_folder, "results/Plots")
    top_poses_folder = os.path.join(output_folder, "results/BestStructs")
    clusters_folder = os.path.join(output_folder, "results/clusters")

    if not os.path.exists(plots_folder):
        os.makedirs(plots_folder)
    if not os.path.exists(top_poses_folder):
        os.makedirs(top_poses_folder)
    if not os.path.exists(clusters_folder):
        os.makedirs(clusters_folder)

    # Plot metrics
    while current_metric <= metrics - 1:
        try:
            analysis.plot_two_metrics(
                total_energy, be, current_metric, output_folder=plots_folder
            )
            analysis.plot_two_metrics(current_metric, be, output_folder=plots_folder)

        except ValueError:
            break
        current_metric += 1

    # Retrieve 100 best structures
    logger.info("Retrieve 100 Best Poses")
    analysis.top_poses(be, 100, top_poses_folder)

    # Clustering of best 2000 best structures
    logger.info(f"Retrieve best cluster poses.")
    if clustering:
        clusters = analysis.cluster_poses(250, be, clusters_folder, nclusts=nclusts)
    if mae:
        sch_python = os.path.join(cs.SCHRODINGER, "utilities/python")
        if not os.path.exists(sch_python):
            sch_python = os.path.join(cs.SCHRODINGER, "run")
        top_poses = glob.glob(os.path.join(top_poses_folder, "*"))
        python_file = os.path.join(cs.DIR, "analysis/to_mae.py")
        for poses in top_poses + clusters:
            command = "{} {} {} --schr {} {}".format(
                sch_python, python_file, poses, cs.SCHRODINGER, "--remove"
            )
            logger.info(command)
            try:
                subprocess.check_call(command.split())
            except ValueError:
                raise ValueError(
                    "Binding energy is not in the default report column (4). Please specify the column by means of "
                    "the next flag. be_column: N "
                )
    plots = glob.glob(os.path.join(plots_folder, "*.png"))
    poses = glob.glob(os.path.join(top_poses_folder, "*"))
    clusters = glob.glob(os.path.join(clusters_folder, "*.png"))
    report = pr.create_report(
        plots,
        clusters,
        poses,
        analysis.best_energies,
        output=os.path.join(output_folder, "summary_results.pdf"),
    )
    return report
