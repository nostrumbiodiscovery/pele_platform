# -*- coding: utf-8 -*-
"""
This is the main module and it is designed to run the PELE platform from
the command-line.
"""

__author__ = "Nostrum Biodiscovery"
__email__ = "pelesupport@nostrumbiodiscovery.com"
__license__ = "Apache-2.0"
__version__ = "1.0"


# Constants
PLOT_TYPE_CHOICES = ('scatter', 'interactive')


# Methods
def parse_args():
    """
    Command line parser.

    Returns
    -------
    csv_path : str
        The path pointing to the input yaml file
    results_folder_path : str
        The path pointing to a PELE Platform's results folder
    plot_type : str
        The plot type to generate
    xdata : list[str]
        The data corresponding to the X axis. It might have a single
        dimension, if only the column number of the metric to plot is
        supplied, or two dimensions if both column number and axis
        label are supplied
    ydata : list[str]
        The data corresponding to the Y axis. It might have a single
        dimension, if only the column number of the metric to plot is
        supplied, or two dimensions if both column number and axis
        label are supplied
    zdata : list[str]
        The data corresponding to the Z axis. It might have a single
        dimension, if only the column number of the metric to plot is
        supplied, or two dimensions if both column number and axis
        label are supplied
    """
    # Parser setup
    from argparse import ArgumentParser

    parser = ArgumentParser(
        description='PELE Platform Plotter')

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-c', '--csv_path', type=str, default=None,
                       help='Path to the input input file')
    group.add_argument('-r', '--results_folder_path', type=str,
                       default=None, help='Path to the input input file')
    parser.add_argument('-t', '--plot_type', type=str,
                        choices=PLOT_TYPE_CHOICES, default=None,
                        help='Plot type. One of ' + str(PLOT_TYPE_CHOICES))
    parser.add_argument("-x", "--xaxis", metavar=('INT', 'STR'),
                        type=str, nargs='*',
                        help="Column number to plot on the X axis. "
                        + "Optionally, the axis label can be defined "
                        + "afterwards", default=None)
    parser.add_argument("-y", "--yaxis",  metavar=('INT', 'STR'),
                        type=str, nargs='*',
                        help="Column number to plot on the Y axis. "
                        + "Optionally, the axis label can be defined "
                        + "afterwards", default=None)
    parser.add_argument("-z", "--zaxis",  metavar=('INT', 'STR'),
                        type=str, nargs='*',
                        help="Column number to plot on the Y axis. "
                        + "Optionally, the axis label can be defined "
                        + "afterwards", default=None)

    # Parse arguments
    parsed_args = parser.parse_args()

    # Extract arguments
    csv_path = parsed_args.csv_path
    results_folder = parsed_args.results_folder_path
    plot_type = parsed_args.plot_type
    xdata = parsed_args.xaxis
    ydata = parsed_args.yaxis
    zdata = parsed_args.zaxis

    return csv_path, results_folder, plot_type, xdata, ydata, zdata


def print_parameters(csv_path, results_folder, plot_type, xdata, ydata, zdata):
    """
    It prints the parameters supplied by the user.

    Parameters
    ----------
    csv_path : str

    """


def parse_axis_data(axis_data):
    """
    It sets the columns and label of the data that wants to be plotted.

    PARAMETERS
    ----------
    axis_data : list[str]
        Axis data to parse

    RETURNS
    -------
    column : tuple[int, str]
        The column number to plot
    label : str
        The label that corresponds to the column metric that will be
        plotted
    """

    if axis_data is None:
        return None, None

    else:
        try:
            column = int(axis_data[0])

            if len(axis_data) > 2:
                raise ValueError
        except ValueError:
            print("Warning: axis data not recognized: {}".format(axis_data))
            return None, None

        if len(axis_data) == 1:
            return column, None

        else:
            metric_name = axis_data[1]

            # Add units to label
            if "energy" in metric_name.lower():
                label = metric_name + " ($kcal/mol$)"
            elif "energies" in metric_name.lower():
                label = metric_name + " ($kcal/mol$)"
            elif "distance" in metric_name.lower():
                label = metric_name + " ($\AA$)"
            elif "rmsd" in metric_name.lower():
                label = metric_name + " ($\AA$)"
            else:
                label = metric_name

            return column, label


# Main workflow to be executed
if __name__ == "__main__":
    # Avoid backend issues in matplotlib
    import matplotlib

    matplotlib.use("Agg")

    # Avoid Python2
    from .Checker.python_version import check_python_version

    check_python_version()

    # Parse command-line arguments
    csv_path, results_folder, plot_type, xdata, ydata, zdata = parse_args()

    # Print header
    from pele_platform.constants import constants
    print(constants.plotter_version_header.format(__version__))

    print_parameters(csv_path, results_folder, plot_type, xdata, ydata, zdata)

    xaxis, xlabel = parse_axis_data(xdata)
    yaxis, ylabel = parse_axis_data(ydata)
    zaxis, zlabel = parse_axis_data(zdata)


def scatter_plot():
    """

    """
    pass


def interactive_plot():
    """

    """
    x_values = []
    y_values = []
    z_values = []
    labels = []
    annotations = []

    with open(reports[0], 'r') as report_file:
        line = report_file.readline()
        if None in x_rows:
            x_rows = [7, ]
            x_name = "RMSD ($\AA$)"
        if None in y_rows:
            y_rows = [5, ]
            y_name = "Energy ($kcal/mol$)"
        if x_name is None:
            x_name = str(line.split("    ")[x_rows[0] - 1])
        if y_name is None:
            y_name = str(line.split("    ")[y_rows[0] - 1])
        if (None not in z_rows) and (z_name is None):
            z_name = str(line.split("    ")[z_rows[0] - 1])
            z_name = addUnits(z_name)

    for report in reports:
        report_directory = os.path.dirname(report)
        report_number = os.path.basename(report).split('_')[-1].split('.')[0]

        with open(report, 'r') as report_file:
            next(report_file)
            for i, line in enumerate(report_file):
                x_total = 0.
                y_total = 0.
                z_total = 0.

                for x_row in x_rows:
                    x_total += float(line.split()[x_row - 1])

                for y_row in y_rows:
                    y_total += float(line.split()[y_row - 1])

                if None not in z_rows:
                    for z_row in z_rows:
                        z_total += float(line.split()[z_row - 1])

                if isnan(x_total) or isnan(y_total) or isnan(z_total):
                    continue

                x_values.append(x_total)
                y_values.append(y_total)
                z_values.append(z_total)

                epoch = report_directory.split('/')[-1]
                if not epoch.isdigit():
                    epoch = '0'

                annotations.append("Epoch: " + epoch + "\n" +
                                   "Trajectory: " + report_number + "\n" +
                                   "Model: " + str(i + 1))

                labels.append(0)

    if z_max is None:
        z_max = max(z_values)

    if z_min is None:
        z_min = min(z_values)

    if z_min == z_max:
        cmap = pyplot.cm.autumn
    else:
        cmap = pyplot.cm.plasma

    norm = pyplot.Normalize(z_min, z_max)

    fig, ax = pyplot.subplots()

    if output_path is not None:
        s = 20
    else:
        s = None

    sc = pyplot.scatter(x_values, y_values, c=z_values, cmap=cmap, s=s,
                        norm=norm)

    ax.margins(0.05)
    ax.set_facecolor('lightgray')
    pyplot.ylabel(y_name)
    pyplot.xlabel(x_name)

    annot = ax.annotate("", xy=(0, 0), xytext=(20, 20),
                        textcoords="offset points",
                        bbox=dict(boxstyle="round", fc="w"),
                        arrowprops=dict(arrowstyle="->"))
    annot.set_visible(False)

    # Activate the colorbar only if the Z axis contains data to plot
    if None not in z_rows:
        cbar = pyplot.colorbar(sc, drawedges=False)
        cbar.ax.set_ylabel(z_name)

    def update_annot(ind):
        """Update the information box of the selected point"""
        pos = sc.get_offsets()[ind["ind"][0]]
        annot.xy = pos
        annot.set_text(annotations[int(ind["ind"][0])])
        annot.get_bbox_patch().set_facecolor(cmap(norm(
            z_values[ind["ind"][0]])))

    def hover(event):
        """Action to perform when hovering the mouse on a point"""
        vis = annot.get_visible()
        if event.inaxes == ax:
            cont, ind = sc.contains(event)
            if cont:
                update_annot(ind)
                annot.set_visible(True)
                fig.canvas.draw_idle()
            else:
                if vis:
                    annot.set_visible(False)
                    fig.canvas.draw_idle()

    # Respond to mouse motion
    fig.canvas.mpl_connect("motion_notify_event", hover)

    # Save or display the plot depending on whether an output path was set or
    # not
    if output_path is not None:
        pyplot.savefig(output_path)
    else:
        pyplot.show()



