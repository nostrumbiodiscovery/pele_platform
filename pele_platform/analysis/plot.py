"""
This module contains classes and methods to construct plots with data
coming from PELE reports.
"""
import os

from pele_platform.Utilities.Helpers.helpers import check_make_folder
from pele_platform.analysis import DataHandler
from pele_platform.constants import constants


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

        check_make_folder(output_folder)

        # Initialize a data handler from the current dataframe and get column names
        metric_to_x, metric_to_y, metric_to_z = self._get_column_names(
            metric_to_x, metric_to_y, metric_to_z)

        # Prepare plot name
        if output_name is None:
            if metric_to_z is not None:
                output_name = "{}_{}_{}_plot.png".format(metric_to_x,
                                                         metric_to_y,
                                                         metric_to_z)
            else:
                output_name = "{}_{}_plot.png".format(metric_to_x, metric_to_y)

        # Replace whitespaces in the output name
        output_name = output_name.replace(" ", "_")
        output_name = os.path.join(output_folder, output_name)

        # Generate plot with matplotlib
        import matplotlib

        matplotlib.use("Agg")
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
            backup_logger(self._logger,
                          "Plotted {} vs {}".format(metric_to_x,
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

        check_make_folder(output_folder)
        metric_to_x, metric_to_y, metric_to_z = \
            self._get_column_names(metric_to_x, metric_to_y)

        # Define output path
        output_name = "{}_{}_kde.png".format(metric_to_x, metric_to_y)
        output_name = output_name.replace(" ", "_")
        output_name = os.path.join(output_folder, output_name)

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
                      clusters, representative_structures=None):
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
        representative_structures : dict[str, tuple[str, int]]
            Dictionary containing the representative structures that
            were selected. Cluster label is the key and value is a list
            with [trajectory, step] of each cluster. If supplied, points
            belonging to representative structures will be represented
        """
        import copy
        from matplotlib.colors import LinearSegmentedColormap

        from pele_platform.analysis.clustering import get_cluster_label
        from pele_platform.Utilities.Helpers.helpers import backup_logger

        check_make_folder(output_folder)

        # Initialize a data handler from the current dataframe
        metric_to_x, metric_to_y, metric_to_z = \
            self._get_column_names(metric_to_x, metric_to_y)

        import matplotlib

        matplotlib.use("Agg")
        from matplotlib import pyplot as plt
        from matplotlib import colors, cm

        # Initialize figure
        fig, ax = plt.subplots(figsize=(6, 6), dpi=100, facecolor="w",
                               edgecolor="k")
        fig.subplots_adjust(right=0.8)  # To make room for the legend

        # Set axis labels
        plt.xlabel(metric_to_x)
        plt.ylabel(metric_to_y)

        # Configurate grid
        ax.set_axisbelow(True)
        ax.grid(True)
        ax.xaxis.grid(color="#AEB6BF", linestyle="dashed")
        ax.yaxis.grid(color="#AEB6BF", linestyle="dashed")
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.spines["left"].set_color("black")
        ax.spines["bottom"].set_color("black")
        ax.set_facecolor("#E6E9EB")

        # Extract the list of cluster labels
        cluster_labels = sorted(list(set(clusters)))

        # Configurate colormap
        if len(cluster_labels) > 18:
            cmap = copy.copy(cm.get_cmap("jet"))
        elif 9 < len(cluster_labels) <= 18:
            cmap = LinearSegmentedColormap.from_list('custom_tab20',
                                                     constants.custom_colors)
        else:
            cmap = copy.copy(cm.get_cmap("Set1"))

        norm = colors.Normalize(vmin=0, vmax=len(cluster_labels))
        cmap.set_under("grey")

        # Values to plot
        all_xs = self._dataframe[metric_to_x]
        all_ys = self._dataframe[metric_to_y]

        colors_used = []
        rep_struct_marker = []
        if representative_structures is not None:
            marker_cm = cm.get_cmap('binary')
            marker_norm = colors.Normalize(vmin=0, vmax=1)
            # Draw points with representative structures
            trajectories = self._dataframe['trajectory']
            steps = self._dataframe['numberOfAcceptedPeleSteps']

            rep_trajs = [traj for (traj, step)
                         in representative_structures.values()]
            rep_steps = [step for (traj, step)
                         in representative_structures.values()]

            for current_cluster in cluster_labels:
                xs = []
                ys = []
                for x, y, cluster, traj, step in zip(all_xs, all_ys,
                                                     clusters,
                                                     trajectories,
                                                     steps):
                    if cluster == current_cluster:
                        traj_idxs = set([i for i, x in enumerate(rep_trajs)
                                         if x == traj])
                        step_idxs = set([i for i, x in enumerate(rep_steps)
                                         if x == step])
                        if len(traj_idxs.intersection(step_idxs)) == 1:
                            sc = ax.scatter([x, ], [y, ], c=[1, ],
                                            zorder=3, marker='x',
                                            cmap=marker_cm,
                                            norm=marker_norm)
                            rep_struct_marker = sc.legend_elements()[0]
                        else:
                            xs.append(x)
                            ys.append(y)

                # In case there is only one point and it is the
                # representative structure
                if len(xs) == 0:
                    continue

                if current_cluster == -1:
                    zorder = 1
                else:
                    zorder = 2

                sc = ax.scatter(xs, ys, c=[current_cluster, ] * len(xs),
                                cmap=cmap, norm=norm, alpha=0.7,
                                zorder=zorder)
                colors_used += sc.legend_elements()[0]

        else:
            # Draw points without representative structures
            for current_cluster in cluster_labels:
                xs = []
                ys = []
                for x, y, cluster in zip(all_xs, all_ys, clusters):
                    if cluster == current_cluster:
                        xs.append(x)
                        ys.append(y)
                if current_cluster == -1:
                    zorder = 1
                else:
                    zorder = 2

                sc = ax.scatter(xs, ys, c=[current_cluster, ] * len(xs),
                                cmap=cmap, norm=norm, alpha=0.7,
                                zorder=zorder)
                colors_used += sc.legend_elements()[0]

        # Configure legend
        cluster_names = []
        for cluster_id in cluster_labels:
            if cluster_id == -1:
                cluster_names.append("Others")
            else:
                cluster_names.append(get_cluster_label(cluster_id))

        if cluster_names[0] == "Others":
            n = cluster_names.pop(0)
            c = colors_used.pop(0)
            cluster_names.append(n)
            colors_used.append(c)

        if len(rep_struct_marker) == 1:
            cluster_names.append("Representative\nstructure")

        ax.legend(colors_used + rep_struct_marker, cluster_names,
                  title="Clusters", loc='center left',
                  bbox_to_anchor=(1, 0.5))

        # Set output name
        if representative_structures is not None:
            output_name = "{}_{}_representatives_plot.png".format(metric_to_x,
                                                                  metric_to_y)
        else:
            output_name = "{}_{}_plot.png".format(metric_to_x, metric_to_y)
        output_name = output_name.replace(" ", "_")
        output_name = os.path.join(output_folder, output_name)

        plt.savefig(output_name, dpi=200, edgecolor="k",
                    orientation="portrait", transparent=True,
                    bbox_inches="tight")

        backup_logger(self._logger,
                      "Plotted {} vs {}".format(metric_to_x, metric_to_y))

    def _get_column_names(self, metric_to_x, metric_to_y, metric_to_z=None):
        """
        Gets column names based on indices.
        Parameters
        ----------
        metric_to_x : int
            Index of column to be plotted on x-axis
        metric_to_y : int
            Index of column to be plotted on y-axis
        metric_to_z : int
            Index of column to be used for the colour bar.
        Returns
        -------
            Column names of x, y and z axes.
        """

        data_handler = DataHandler.from_dataframe(self._dataframe)
        # Ensure that metrics are strings pointing to dataframe columns
        if str(metric_to_x).isdigit():
            metric_to_x = data_handler.get_column_name(metric_to_x)
        if str(metric_to_y).isdigit():
            metric_to_y = data_handler.get_column_name(metric_to_y)
        if metric_to_z is not None and str(metric_to_z).isdigit():
            metric_to_z = data_handler.get_column_name(metric_to_z)
        return metric_to_x, metric_to_y, metric_to_z
