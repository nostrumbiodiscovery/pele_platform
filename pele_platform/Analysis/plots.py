import os
import pandas as pd
import matplotlib.pyplot as plt
import AdaptivePELE.analysis.selectOnPlot as sp
import AdaptivePELE.analysis.bestStructs as bs

EPOCH = "epoch"
STEPS = 3
TRAJECTORY = "trajectory"



class PostProcessor():

    def __init__(self, report_name, traj_name, simulation_path, topology=False):
        self.report_name = report_name
        self.traj_name = traj_name
        self.simulation_path = simulation_path
        self.data = self.retrive_data()
        self.topology = False

    def retrive_data(self, separator=","):
        summary_csv_filename = os.path.join(self.simulation_path, "summary.csv")
        if not os.path.exists(summary_csv_filename):
            sp.concat_reports_in_csv(adaptive_results_path=self.simulation_path, output_file_path=summary_csv_filename,
                                  report_prefix=self.report_name, trajectory_prefix=self.traj_name, separator_out=separator)
        dataframe = pd.read_csv(summary_csv_filename, sep=separator, engine='python', header=0)
        return dataframe

    def plot_two_metrics(self, column_to_x, column_to_y, column_to_z=False, output_name=None):

        column_to_x = column_to_x if not str(column_to_x).isdigit() else self._get_column_name(self.data, column_to_x)
        column_to_y = column_to_y if not str(column_to_y).isdigit() else self._get_column_name(self.data, column_to_y)

        fig, ax = plt.subplots()
        if column_to_z:
            column_to_z = column_to_z if not str(column_to_z).isdigit() else self._get_column_name(self.data, column_to_z)
            output_name = output_name if output_name else "{}_{}_{}_plot.png".format(column_to_x, column_to_y, column_to_z)
            output_name = output_name.replace(" ", "_")
            pts = ax.scatter(self.data[column_to_x], self.data[column_to_y], c=self.data[column_to_z], s=20)
            cbar = plt.colorbar(pts)
            cbar.ax.set_ylabel(column_to_z)
            ax.set_xlabel(column_to_x)
            ax.set_ylabel(column_to_y)
            plt.savefig(output_name)
            print("Plotted {} vs {} vs {}".format(column_to_x, column_to_y, column_to_z))
        else:
            output_name = output_name if output_name else "{}_{}_plot.png".format(column_to_x, column_to_y)
            output_name = output_name.replace(" ", "_")
            pts = ax.scatter(self.data[column_to_x], self.data[column_to_y], s=20)
            plt.savefig(output_name)
            print("Plotted {} vs {}".format(column_to_x, column_to_y))

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

    def _get_column_name(self, df, column_digit):
        return list(df)[int(column_digit)-1]


def analyse_simulation(report_name, traj_name, simulation_path):
    analysis = PostProcessor(report_name, traj_name, simulation_path)
    metrics = len(list(analysis.data)) - 1 #Discard epoch as metric
    sasa = 6
    be = 5
    total_energy = 4
    current_metric = sasa
    while current_metric <= metrics-1:
        try:
            analysis.plot_two_metrics(be, total_energy, current_metric)
            analysis.plot_two_metrics(be, current_metric)
        except ValueError:
            break
        current_metric += 1

    analysis.top_poses(be, 100)


if __name__ == "__main__":
    analyse_simulation("*report*", "*trajectory*", "/work/NBD_Utilities/PELE/PELE_Softwares/pele_platform/tests/STR_Pele/output/")


