import os
import glob
import subprocess
import mdtraj
import numpy as np
import pandas as pd
from sklearn import mixture
from multiprocessing import Pool
import matplotlib.pyplot as plt
import AdaptivePELE.analysis.selectOnPlot as sp
import AdaptivePELE.analysis.bestStructs as bs
from pele_platform.constants import constants as cs

EPOCH = "epoch"
STEPS = 3
TRAJECTORY = "trajectory"


def _extract_coords(info):
    p, v, resname = info
    # Most time consuming step 0.1
    traj = mdtraj.load_frame(p, v, top="topology.pdb")
    atoms_info = traj.topology.to_dataframe()[0]
    condition = atoms_info['resName'] == resname
    atom_numbers_ligand = atoms_info[condition].index.tolist()
    coords = []
    for atom_num in atom_numbers_ligand:
        try:
            coords.extend(traj.xyz[0][atom_num].tolist())
        except IndexError:
            continue
    return np.array(coords).ravel()


class PostProcessor():

    def __init__(self, report_name, traj_name, simulation_path, cpus, topology=False, residue=False):
        self.report_name = report_name
        self.traj_name = traj_name
        self.simulation_path = simulation_path
        self.data = self.retrive_data()
        self.topology = topology
        self.residue = residue
        self.cpus = cpus

    def retrive_data(self, separator=","):
        summary_csv_filename = os.path.join(self.simulation_path, "summary.csv")
        if not os.path.exists(summary_csv_filename):
            sp.concat_reports_in_csv(adaptive_results_path=self.simulation_path, output_file_path=summary_csv_filename,
                                  report_prefix=self.report_name, trajectory_prefix=self.traj_name, separator_out=separator)
        dataframe = pd.read_csv(summary_csv_filename, sep=separator, engine='python', header=0)
        return dataframe

    def plot_two_metrics(self, column_to_x, column_to_y, column_to_z=False, output_name=None, output_folder="."):
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

        column_to_x = column_to_x if not str(column_to_x).isdigit() else self._get_column_name(self.data, column_to_x)
        column_to_y = column_to_y if not str(column_to_y).isdigit() else self._get_column_name(self.data, column_to_y)

        fig, ax = plt.subplots()
        if column_to_z:
            column_to_z = column_to_z if not str(column_to_z).isdigit() else self._get_column_name(self.data, column_to_z)
            output_name = output_name if output_name else "{}_{}_{}_plot.png".format(column_to_x, column_to_y, column_to_z)
            output_name = os.path.join(output_folder,output_name).replace(" ", "_")
            pts = ax.scatter(self.data[column_to_x], self.data[column_to_y], c=self.data[column_to_z], s=20)
            cbar = plt.colorbar(pts)
            cbar.ax.set_ylabel(column_to_z)
            ax.set_xlabel(column_to_x)
            ax.set_ylabel(column_to_y)
            plt.savefig(output_name)
            print("Plotted {} vs {} vs {}".format(column_to_x, column_to_y, column_to_z))
        else:
            output_name = output_name if output_name else "{}_{}_plot.png".format(column_to_x, column_to_y)
            output_name = os.path.join(output_folder,output_name).replace(" ", "_")
            pts = ax.scatter(self.data[column_to_x], self.data[column_to_y], s=20)
            ax.set_xlabel(column_to_x)
            ax.set_ylabel(column_to_y)
            plt.savefig(output_name)
            print("Plotted {} vs {}".format(column_to_x, column_to_y))
        return output_name

    def top_poses(self, metric, n_structs, output="BestStructs"):
        metric = metric if not str(metric).isdigit() else self._get_column_name(self.data, metric)
        best_poses = self.data.nsmallest(n_structs, metric)
        self._extract_poses(best_poses, metric, output)

    def _extract_poses(self, poses, metric, output):
        #Process data to be outputted
        values = poses[metric].tolist()
        paths = poses[TRAJECTORY].tolist()
        epochs = poses[EPOCH].tolist()
        file_ids = [traj[:-4].split("_")[-1] for traj in paths]
        steps = self._get_column_name(self.data, STEPS)
        step_indexes = poses[steps].tolist()
        out_freq = 1
        files_out = ["epoch{}_trajectory_{}.{}_{}{}.pdb".format(epoch, report, int(step), metric.replace(" ",""), value) \
           for epoch, step, report, value in zip(epochs, step_indexes, file_ids, values)]

        #Read traj and output sanpshot
        for f_id, f_out, step, path in zip(file_ids, files_out, step_indexes, paths):
            if not self.topology:
                bs.extract_snapshot_from_pdb(path, f_id, output, self.topology, step, out_freq, f_out)
            else:
                bs.extract_snapshot_from_xtc(path, f_id, output, self.topology, step, out_freq, f_out)

    def cluster_poses(self, n_structs, metric, output):
        assert self.residue, "Set residue ligand name to clusterize"
        metric = metric if not str(metric).isdigit() else self._get_column_name(self.data, metric)
        best_poses = self.data.nsmallest(n_structs, metric)
        self._cluster(best_poses, metric, output)

    def _cluster(self, poses, metric, output):
        # Extract metric values
        values = poses[metric].tolist()
        epochs = poses[EPOCH].tolist()
        paths = poses[TRAJECTORY].tolist()
        file_ids = [traj[:-4].split("_")[-1] for traj in paths]
        steps = self._get_column_name(self.data, STEPS)
        snapshots = poses[steps].tolist()
        # Extract coords 
        all_coords = []
        pool = Pool(processes=self.cpus)
        input_pool = [[p,v,self.residue] for p, v in zip(paths, snapshots)]
        all_coords = pool.map(_extract_coords, input_pool)
        # Clusterize
        assert all_coords[0][0], "Ligand not found check the option --resname. i.e python interactive.py 5 6 7 --resname LIG"
        try:
            clf = mixture.GaussianMixture(n_components=10, covariance_type='full')
            labels = clf.fit_predict(all_coords)
            indexes = labels
        except ValueError:
            indexes = [1]
        n_clusters = len(set(indexes))
        files_out = ["cluster{}_epoch{}_trajectory_{}.{}_{}{}.pdb".format(cluster, epoch, report, int(step), metric.replace(" ",""), value) \
           for epoch, step, report, value, cluster in zip(epochs, snapshots, file_ids, values, indexes)]
        all_metrics = []
        for n_cluster in range(n_clusters-1):
            metrics = {value:idx for idx, (value, cluster) in enumerate(zip(values, indexes)) if n_cluster == cluster}
            out_freq = 1
            cluster_metrics = list(metrics.keys())
            max_idx = metrics[np.min(cluster_metrics)]
            max_traj = paths[max_idx]
            max_snapshot = snapshots[max_idx]
            output_traj = files_out[max_idx]
            input_traj = file_ids[max_idx]
            if not self.topology:
                bs.extract_snapshot_from_pdb(max_traj, input_traj, output, self.topology, max_snapshot, out_freq, output_traj)
            else:
                bs.extract_snapshot_from_xtc(max_traj, input_traj, output, self.topology, max_snapshot, out_freq, output_traj)
            all_metrics.append(cluster_metrics)
        fig, ax = plt.subplots()
        try:
            ax.boxplot(all_metrics)
        except IndexError:
            print("Samples to disperse to produce a cluster")
            return
        ax.set_ylabel(metric)
        ax.set_xlabel("Cluster number")
        plt.savefig(os.path.join(output, "clusters_{}_boxplot.png".format(metric)))




    def _get_column_name(self, df, column_digit):
        return list(df)[int(column_digit)-1]


def analyse_simulation(report_name, traj_name, simulation_path, residue, output_folder=".", cpus=5, clustering=True, mae=False):
    analysis = PostProcessor(report_name, traj_name, simulation_path, cpus, residue=residue)

    metrics = len(list(analysis.data)) - 1 #Discard epoch as metric
    sasa = 6
    be = 5
    total_energy = 4
    current_metric = sasa
    plots_folder = os.path.join(output_folder, "results/Plots")
    top_poses_folder = os.path.join(output_folder, "results/BestStructs")
    clusters_folder = os.path.join(output_folder, "results/clusters")

    os.makedirs(plots_folder, exist_ok=True)
    os.makedirs(top_poses_folder, exist_ok=True)
    os.makedirs(clusters_folder, exist_ok=True)


    # Plot metrics
    while current_metric <= metrics-1:
        try:
            analysis.plot_two_metrics(total_energy, be, current_metric, output_folder=plots_folder)
            analysis.plot_two_metrics(current_metric, be, output_folder=plots_folder)
        except ValueError:
            break
        current_metric += 1

    #Retrieve 100 best structures
    print("Retrieve 100 Best Poses")
    analysis.top_poses(be, 100, top_poses_folder)

    #Clustering of best 2000 best structures
    print("Retrieve 10 best cluster poses")
    if clustering:
        analysis.cluster_poses(250, be, clusters_folder)

    if mae:
        sch_python = os.path.join(cs.SCHRODINGER, "utilities/python")
        if not os.path.exists(sch_python):
            sch_python = os.path.join(cs.SCHRODINGER, "run")
        top_poses = glob.glob(os.path.join(top_poses_folder, "*"))
        python_file = os.path.join(cs.DIR, "Analysis/to_mae.py")
        for poses in top_poses:
            command = "{} {} {} --schr {} {}".format(sch_python, python_file, poses,  cs.SCHRODINGER, "--remove") 
            print(command)
            subprocess.call(command.split())
    plots = glob.glob(os.path.join(plots_folder, "*.png"))
    poses = glob.glob(os.path.join(top_poses_folder, "*"))
    clusters = glob.glob(os.path.join(clusters_folder, "*.png"))
    return plots, poses, clusters
       




if __name__ == "__main__":
    analyse_simulation("report_", "trajectory_", "/scratch/jobs/dsoler/water/trypsin/BEN_Pele_36/output/", residue="BEN")
