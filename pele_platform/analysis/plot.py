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
                       self._dataframe[metric_to_y])
            ax.set_xlabel(metric_to_x)
            ax.set_ylabel(metric_to_y)
            plt.savefig(output_name)
            backup_logger(self._logger, "Plotted {} vs {}".format(metric_to_x,
                                                                  metric_to_y))
        return output_name

    def plot_kde(self, metric_to_x, metric_to_y, output_folder, kde_structs):
        """
        Given 2 metrics, it generates the kde plot.

        Parameters
        ----------
        metric_to_x : str or int
            The metric id to plot in the X axis. It can be either a string
            with the name of the metric or an integer with the column index
        metric_to_y : str or int
            The metric id to plot in the Y axis. It can be either a string
            with the name of the metric or an integer with the column index
        output_folder : str
            The path where the plot will be saved
        kde_structs : int
            The number of structures to represent in the plot
        """
        import seaborn as sb

        data_handler = DataHandler.from_dataframe(self._dataframe)

        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

        # Ensure that metrics are strings pointing to dataframe columns
        if str(metric_to_x).isdigit():
            metric_to_x = data_handler._get_column_name(metric_to_x)
        if str(metric_to_y).isdigit():
            metric_to_y = data_handler._get_column_name(metric_to_y)

        # Define output path
        output_name = "{}_{}_kde.png".format(metric_to_x,metric_to_y)
        output_name = os.path.join(output_folder, output_name)
        output_name = output_name.replace(" ", "_")

        # Filter out the number of structures from dataframe to plot
        structures_to_keep = min(int(kde_structs), len(self._dataframe) - 1)
        sorted_df = self._dataframe.sort_values(metric_to_y, ascending=True)
        top = sorted_df[0:structures_to_keep]

        # Plot and save it
        plot = sb.kdeplot(top[metric_to_x], top[metric_to_y],
                          cmap="crest", fill=False,
                          shade=True, cbar=True)
        figure = plot.get_figure()
        figure.savefig(output_name)
        return output_name

    def plot_clusters(self, metric_to_x, metric_to_y, output_folder,
                      clusters):
        """
        It creates a scatter plot with the two metrics that are supplied
        and displays the points belonging to each top cluster with a
        different color.

        Parameters
        ----------
        metric_to_x : str or int
            The metric id to plot in the X axis. It can be either a string
            with the name of the metric or an integer with the column index
        metric_to_y : str or int
            The metric id to plot in the Y axis. It can be either a string
            with the name of the metric or an integer with the column index
        output_folder : str
            The path where the plot will be saved
        clusters : a numpy.array object
            The array of cluster labels that were obtained
        """
        import copy
        from string import ascii_uppercase
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

        import matplotlib
        matplotlib.use('Agg')
        from matplotlib import pyplot as plt
        from matplotlib import colors, cm

        # Initialize figure
        fig, ax = plt.subplots(figsize=(6, 6), dpi=100, facecolor='w',
                               edgecolor='k')

        # Set axis labels
        plt.xlabel(metric_to_x)
        plt.ylabel(metric_to_y)

        # Configurate grid
        ax.set_axisbelow(True)
        ax.grid(True)
        ax.xaxis.grid(color='#AEB6BF', linestyle='dashed')
        ax.yaxis.grid(color='#AEB6BF', linestyle='dashed')
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_color('black')
        ax.spines['bottom'].set_color('black')
        ax.set_facecolor('#E6E9EB')

        # Extract the list of cluster labels
        cluster_labels = sorted(list(set(clusters)))

        # Configurate colormap
        cmap = copy.copy(cm.get_cmap("Set1"))
        norm = colors.Normalize(vmin=0, vmax=len(cluster_labels))
        cmap.set_under('grey')

        # Values to plot
        all_xs = self._dataframe[metric_to_x]
        all_ys = self._dataframe[metric_to_y]

        # Draw points
        colors_used = []
        for current_cluster in cluster_labels:
            xs = []
            ys = []
            for x, y, cluster in zip(all_xs, all_ys, clusters):
                if (cluster == current_cluster):
                    xs.append(x)
                    ys.append(y)
            if (current_cluster == -1):
                zorder = 1
            else:
                zorder = 2
            sc = ax.scatter(xs, ys, c=[current_cluster, ] * len(xs), cmap=cmap,
                            norm=norm, alpha=0.7, zorder=zorder)
            colors_used += sc.legend_elements()[0]

        # Configurate legend
        cluster_names = []
        for cluster_id in cluster_labels:
            if cluster_id == -1:
                cluster_names.append('Others')
            else:
                cluster_names.append(ascii_uppercase[cluster_id])

        if cluster_names[0] == 'Others':
            n = cluster_names.pop(0)
            c = colors_used.pop(0)
            cluster_names.append(n)
            colors_used.append(c)
        legend1 = ax.legend(colors_used, cluster_names, title="Clusters")

        # Set output name
        output_name = "{}_{}_plot.png".format(metric_to_x, metric_to_y)
        output_name = os.path.join(output_folder, output_name).replace(" ", "_")

        plt.savefig(output_name, dpi=200, edgecolor='k',
                    orientation='portrait', transparent=True)

        backup_logger(self._logger, "Plotted {} vs {}".format(metric_to_x,
                                                              metric_to_y))
