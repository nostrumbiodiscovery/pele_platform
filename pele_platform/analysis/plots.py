import os
import shutil
import glob
import subprocess
import mdtraj
import numpy as np
import pandas as pd
from sklearn import mixture, cluster
from multiprocessing import Pool
import matplotlib.pyplot as plt
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


def _extract_raw_coords(info):
    """
    This method attempts to be an alternative coordinate extraction with
    fewer memory demands.

    .. todo ::
       * info parameter should be replaced by something more self-explanatory.

    Parameters
    ----------
    info : tuple[str, int, str, str]
        The information
    """
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


class PostProcessor:
    def __init__(
            self,
            report_name,
            traj_name,
            simulation_path,
            cpus,
            topology=False,
            residue=False,
            be_column=4,
            limit_column=6,
            te_column=3,
            clustering_method="GaussianMixture",
            bandwidth=None,
            logger=None,
    ):
        self.report_name = report_name
        self.traj_name = traj_name
        self.simulation_path = simulation_path
        self.data = self.retrive_data()
        self.topology = topology
        self.residue = residue
        self.be_column = be_column
        self.limit_column = limit_column
        self.te_column = te_column
        self.cpus = cpus
        self.clustering_method = clustering_method
        self.bandwidth = bandwidth
        self.logger = logger

    def retrive_data(self, separator=","):
        summary_csv_filename = os.path.join(self.simulation_path, "summary.csv")
        if not os.path.exists(summary_csv_filename) or self._moved_folder(
                summary_csv_filename
        ):
            try:
                sp.concat_reports_in_csv(
                    adaptive_results_path=self.simulation_path,
                    output_file_path=summary_csv_filename,
                    report_prefix=self.report_name,
                    trajectory_prefix=self.traj_name,
                    separator_out=separator,
                )
            except ValueError:
                raise ValueError(
                    "Not report found under {}. Did you point to the right folder?".format(
                        os.getcwd()
                    )
                )
        dataframe = pd.read_csv(
            summary_csv_filename, sep=separator, header=0, float_precision="round_trip"
        )
        dataframe_filtered = self._remove_outliers(dataframe)
        return dataframe_filtered

    def _remove_outliers(self, dataframe):

        cols = list(dataframe.columns)
        n_points_to_remove = int(
            len(dataframe[cols[0]]) * 0.02
        )  # remove top 2% of each energy
        dataframe_filtered = dataframe.sort_values(cols[3], ascending=False).iloc[
                             n_points_to_remove:
                             ]
        dataframe_filtered = dataframe_filtered.sort_values(
            cols[4], ascending=False
        ).iloc[n_points_to_remove:]

        return dataframe_filtered

    def plot_two_metrics(
            self,
            column_to_x,
            column_to_y,
            column_to_z=None,
            output_name=None,
            output_folder=".",
    ):
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

        column_to_x = (
            column_to_x
            if not str(column_to_x).isdigit()
            else self._get_column_name(self.data, column_to_x)
        )
        column_to_y = (
            column_to_y
            if not str(column_to_y).isdigit()
            else self._get_column_name(self.data, column_to_y)
        )

        fig, ax = plt.subplots()
        if column_to_z:
            column_to_z = (
                column_to_z
                if not str(column_to_z).isdigit()
                else self._get_column_name(self.data, column_to_z)
            )
            output_name = (
                output_name
                if output_name
                else "{}_{}_{}_plot.png".format(column_to_x, column_to_y, column_to_z)
            )
            output_name = os.path.join(output_folder, output_name).replace(" ", "_")

            pts = ax.scatter(
                self.data[column_to_x],
                self.data[column_to_y],
                c=self.data[column_to_z],
                s=20,
            )
            cbar = plt.colorbar(pts)
            cbar.ax.set_ylabel(column_to_z)
            ax.set_xlabel(column_to_x)
            ax.set_ylabel(column_to_y)
            plt.savefig(output_name)
            backup_logger(
                self.logger,
                "Plotted {} vs {} vs {}".format(column_to_x, column_to_y, column_to_z),
            )
        else:
            output_name = (
                output_name
                if output_name
                else "{}_{}_plot.png".format(column_to_x, column_to_y)
            )
            output_name = os.path.join(output_folder, output_name).replace(" ", "_")
            ax.scatter(self.data[column_to_x], self.data[column_to_y], s=20)
            ax.set_xlabel(column_to_x)
            ax.set_ylabel(column_to_y)
            plt.savefig(output_name)
            backup_logger(
                self.logger, "Plotted {} vs {}".format(column_to_x, column_to_y)
            )
        return output_name

    def top_poses(self, metric, n_structs, output="BestStructs"):
        metric = (
            metric
            if not str(metric).isdigit()
            else self._get_column_name(self.data, metric)
        )
        best_poses = self.data.nsmallest(n_structs, metric)
        self._extract_poses(best_poses, metric, output)

    def _extract_poses(self, poses, metric, output):
        # Process data to be outputted
        values = poses[metric].tolist()
        self.best_energies = values
        paths = poses[TRAJECTORY].tolist()
        epochs = poses[EPOCH].tolist()
        file_ids = [traj[:-4].split("_")[-1] for traj in paths]
        steps = self._get_column_name(self.data, STEPS)
        step_indexes = poses[steps].tolist()
        out_freq = 1

        distance_key = "distance0.5"
        if distance_key in poses.columns:
            dist_values = poses[distance_key].tolist()
            filename_template = "{}.{}.{}_BindEner{:.2f}_AtomDist{:.2f}.pdb"
            files_out = [
                filename_template.format(epoch, report, int(step), value, dist)
                for epoch, step, report, value, dist in zip(
                    epochs, step_indexes, file_ids, values, dist_values
                )
            ]
        else:
            filename_template = "{}.{}.{}_BindEner{:.2f}.pdb"
            files_out = [
                filename_template.format(epoch, report, int(step), value)
                for epoch, step, report, value in zip(
                    epochs, step_indexes, file_ids, values
                )
            ]

        # Read trajectory and output snapshot
        for f_id, f_out, step, path in zip(file_ids, files_out, step_indexes, paths):
            if not self.topology:
                try:
                    bs.extract_snapshot_from_pdb(
                        path,
                        f_id,
                        output,
                        self.topology,
                        step,
                        out_freq,
                        f_out,
                        self.logger,
                    )
                except UnicodeDecodeError:
                    raise Exception(
                        "XTC output being treated as PDB. Please specify XTC with the next flag. traj: "
                        "'trajectory_name.xtc' in the input.yaml "
                    )
            else:
                bs.extract_snapshot_from_xtc(
                    path,
                    f_id,
                    output,
                    self.topology,
                    step,
                    out_freq,
                    f_out,
                    self.logger,
                )

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
