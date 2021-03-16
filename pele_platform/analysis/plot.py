"""
This module contains classes and methods to construct plots with data
coming from PELE reports.
"""
import os
from pele_platform.analysis import DataHandler


# TODO this method should be moved to another module
def _extract_coords(info):
    import numpy as np
    import mdtraj

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
                         output_name=None, output_folder=".", colors=None, limit_column=6):
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
        colors : list
            List of cluster indices.
        """
        from pele_platform.Utilities.Helpers.helpers import backup_logger

        # Initialize a data handler from the current dataframe
        data_handler = DataHandler.from_dataframe(self._dataframe)

        # Ensure that output_folder exists
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

        # Ensure that metrics are strings pointing to dataframe columns
        if str(metric_to_x).isdigit():
            metric_to_x = data_handler._get_column_name(metric_to_x)
        if str(metric_to_y).isdigit():
            metric_to_y = data_handler._get_column_name(metric_to_y)
        if metric_to_z is not None and str(metric_to_z).isdigit():
            metric_to_z = data_handler._get_column_name(metric_to_z)

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
        import matplotlib
        matplotlib.use('Agg')
        from matplotlib import pyplot as plt

        fig, ax = plt.subplots()
        if metric_to_z is not None:
            colors = colors if colors is not None else self._dataframe[metric_to_z]
            scatter = ax.scatter(self._dataframe[metric_to_x],
                                 self._dataframe[metric_to_y],
                                 c=colors,
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
            colors = colors if colors is not None else None
            ax.scatter(self._dataframe[metric_to_x],
                       self._dataframe[metric_to_y],
                       c=colors,
                       s=20)
            ax.set_xlabel(metric_to_x)
            ax.set_ylabel(metric_to_y)
            plt.savefig(output_name)
            backup_logger(self._logger, "Plotted {} vs {}".format(metric_to_x,
                                                                  metric_to_y))
        return output_name

    def plot_kde(self, column_to_x, column_to_y, output_folder, kde_structs):
        import seaborn as sb

        data_handler = DataHandler.from_dataframe(self._dataframe)

        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

        column_to_x = column_to_x if not str(column_to_x).isdigit() else data_handler._get_column_name(column_to_x)
        column_to_y = column_to_y if not str(column_to_y).isdigit() else data_handler._get_column_name(column_to_y)

        output_name = "{}_{}_kde.png".format(column_to_x, column_to_y)
        output_name = os.path.join(output_folder,output_name).replace(" ", "_")
        top_1000 = self._dataframe.sort_values(column_to_y, ascending=True)[0:int(kde_structs)]
        plot = sb.kdeplot(top_1000[column_to_x], top_1000[column_to_y], cmap="crest", fill=False, shade=True, cbar=True)
        figure = plot.get_figure()
        figure.savefig(output_name)

        return output_name
