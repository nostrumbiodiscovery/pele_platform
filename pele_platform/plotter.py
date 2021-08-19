# -*- coding: utf-8 -*-
"""
This is the main module and it is designed to run the PELE platform from
the command-line.
"""

__author__ = "Nostrum Biodiscovery"
__email__ = "pelesupport@nostrumbiodiscovery.com"
__license__ = "Apache-2.0"
__version__ = "1.0"


# General imports
import os


# Constants
PLOT_TYPE_CHOICES = ('scatter', 'interactive')
DEFAULT_PLOT_TYPE = 'scatter'


# Classes
class PlotData(object):
    """
    It represents the data of all the axes of a plot.
    """

    def __init__(self, x_data, y_data, z_data=None):
        """
        It initializes a PlotData object.

        Parameters
        ----------
        x_data : an AxisData object
            The AxisData object that corresponds to the X axis
        y_data : an AxisData object
            The AxisData object that corresponds to the Y axis
        z_data : an AxisData object
            The AxisData object that corresponds to the Z axis. It
            is optional
        """
        self.x_data = x_data
        self.y_data = y_data
        self.z_data = z_data

    def get_xs_from_pele_data(self, pele_data):
        """
        It returns the X values to plot from a DataFrame containing
        PELE data.

        Parameters
        ----------
        pele_data : a Pandas.DataFrame object
            The DataFrame containing the PELE simulation data from
            where the values will be extracted
        """
        if self.x_data.category == 'report_column':
            return pele_data.iloc[:, self.x_data.column - 1]
        else:
            raise NotImplementedError()

    def get_ys_from_pele_data(self, pele_data):
        """
        It returns the Y values to plot from a DataFrame containing
        PELE data.

        Parameters
        ----------
        pele_data : a Pandas.DataFrame object
            The DataFrame containing the PELE simulation data from
            where the values will be extracted
        """
        if self.y_data.category == 'report_column':
            return pele_data.iloc[:, self.y_data.column - 1]
        else:
            raise NotImplementedError()

    def get_zs_from_pele_data(self, pele_data):
        """
        It returns the Z values to plot from a DataFrame containing
        PELE data.

        Parameters
        ----------
        pele_data : a Pandas.DataFrame object
            The DataFrame containing the PELE simulation data from
            where the values will be extracted
        """
        if isinstance(self.z_data, EmptyAxisData):
            return None

        if self.z_data.category == 'report_column':
            return pele_data.iloc[:, self.z_data.column - 1]
        else:
            raise NotImplementedError()


class AxisData(object):
    """
    It represents the data to plot in an axis.
    """

    _CATEGORIES = ['report_column', ]

    def __init__(self, label, column=None, category='report_column',
                 lowest=None, highest=None):
        """
        It initializes an AxisData object.

        Parameters
        ----------
        label : str
            The label corresponding to the current axis data
        column : int
            If metric is found in PELE's report files, it specifies
            the column that corresponds to the metric to plot in
            the current axis, starting at index 1
        category : str
            The type of metric that will be plotted in the current axis.
            Default is `report_column` which implies that data will be
            obtained from PELE's report files. Currently, this is the
            only data category that is supported
        lowest : float
            The lowest value to plot in the current axis. Default is
            'None', which does not apply any constraint
        highest : float
            The highest value to plot in the current axis. Default is
            `None`, which does not apply any constraint
        """
        self.label = label
        self.column = column

        if category not in self._CATEGORIES:
            raise ValueError('Wrong AxisData category:', category)

        self.category = category
        self.lowest = lowest
        self.highest = highest


class EmptyAxisData(AxisData):
    """
    It represents an empty AxisData.
    """

    def __init__(self):
        """
        It initializes an EmptyAxisData object.
        """
        super().__init__(None)


class UnlabelledAxisData(AxisData):
    """
    It represents an unlabelled AxisData.
    """

    def __init__(self, column, category='report_column',
                 lowest=None, highest=None):
        """
        It initializes an EmptyAxisData object.

        Parameters
        ----------
        column : int
            If metric is found in PELE's report files, it specifies
            the column that corresponds to the metric to plot in
            the current axis, starting at index 1
        category : str
            The type of metric that will be plotted in the current axis.
            Default is `report_column` which implies that data will be
            obtained from PELE's report files. Currently, this is the
            only data category that is supported
        lowest : float
            The lowest value to plot in the current axis. Default is
            'None', which does not apply any constraint
        highest : float
            The highest value to plot in the current axis. Default is
            `None`, which does not apply any constraint
        """
        super().__init__(None, column=column, category=category,
                         lowest=lowest, highest=highest)


# Methods
def parse_args():
    """
    Command line parser.

    Returns
    -------
    csv_file : str
        The path pointing to the PELE Platform's csv file
    results_folder : str
        The path pointing to a PELE Platform's results folder
    output_folder : str
        The path pointing to a PELE Platform's output folder
    report_name : str
        The name of PELE's report files. Only required when extracting
        data from PELE's output folder
    trajectory_name : str
        The name of PELE's trajectory files. Only required when
        extracting data from PELE's output folder
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
    group.add_argument('-c', '--csv_file', type=str, default=None,
                       help='Path to PELE csv file')
    group.add_argument('-r', '--results_folder', type=str,
                       default=None, help='Path to results folder')
    group.add_argument('-o', '--output_folder', type=str,
                       default=None, help='Path to PELE output folder')
    parser.add_argument('--report_name', type=str, default='report',
                        help='Name of PELE\'s report files. ' +
                        'Only required when extracting data from ' +
                        'PELE\'s output folder')
    parser.add_argument('--trajectory_name', type=str,
                        default='trajectory.pdb',
                        help='Name of PELE\'s trajectory files. ' +
                        'Only required when extracting data from ' +
                        'PELE\'s output folder')
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
    parser.add_argument('--xlowest', type=float, default=None,
                        help='Lowest value to plot on the X axis')
    parser.add_argument('--xhighest', type=float, default=None,
                        help='Highest value to plot on the X axis')
    parser.add_argument('--ylowest', type=float, default=None,
                        help='Lowest value to plot on the Y axis')
    parser.add_argument('--yhighest', type=float, default=None,
                        help='Highest value to plot on the Y axis')
    parser.add_argument('--zlowest', type=float, default=None,
                        help='Lowest value to plot on the Z axis')
    parser.add_argument('--zhighest', type=float, default=None,
                        help='Highest value to plot on the Z axis')

    # Parse arguments
    parsed_args = parser.parse_args()

    # Extract arguments
    csv_file = parsed_args.csv_file
    results_folder = parsed_args.results_folder
    output_folder = parsed_args.output_folder
    report_name = parsed_args.report_name
    trajectory_name = parsed_args.trajectory_name
    plot_type = parsed_args.plot_type
    xdata = parsed_args.xaxis
    ydata = parsed_args.yaxis
    zdata = parsed_args.zaxis
    xlowest = parsed_args.xlowest
    xhighest = parsed_args.xhighest
    ylowest = parsed_args.ylowest
    yhighest = parsed_args.yhighest
    zlowest = parsed_args.zlowest
    zhighest = parsed_args.zhighest

    # Check parameters
    if csv_file is not None:
        if not os.path.isfile(csv_file):
            raise ValueError('CSV file not found at:', csv_file)

    if results_folder is not None:
        if (not os.path.isdir(results_folder) and
                os.path.isfile(os.path.join(results_folder, 'data.csv'))):
            raise ValueError('Wrong PELE results folder at:', results_folder)

    if output_folder is not None:
        if not os.path.isdir(output_folder):
            raise ValueError('Wrong PELE output folder at:', output_folder)

    return csv_file, results_folder, output_folder, \
        report_name, trajectory_name, plot_type, xdata, ydata, zdata, \
        xlowest, xhighest, ylowest, yhighest, zlowest, zhighest


def print_parameters(csv_file, results_folder, output_folder,
                     report_name, trajectory_name,
                     plot_type, xdata, ydata, zdata):
    """
    It prints the parameters supplied by the user.

    Parameters
    ----------
    csv_file : str
        The path pointing to the PELE Platform's csv file
    results_folder : str
        The path pointing to a PELE Platform's results folder
    output_folder : str
        The path pointing to a PELE Platform's output folder
    report_name : str
        The name of PELE's report files. Only required when extracting
        data from PELE's output folder
    trajectory_name : str
        The name of PELE's trajectory files. Only required when
        extracting data from PELE's output folder
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

    print('Input parameters:')
    print(' - csv_file:', csv_file)
    print(' - results_folder:', results_folder)
    print(' - output_folder:', output_folder)
    print(' - report_name:', report_name)
    print(' - trajectory_name:', trajectory_name)
    print(' - plot_type:', plot_type)
    print(' - xdata:', xdata)
    print(' - ydata:', ydata)
    print(' - zdata:', zdata)


def parse_axis_data(axis_data, lowest, highest):
    """
    It sets the columns and label of the data that wants to be plotted.

    PARAMETERS
    ----------
    axis_data : list[str]
        Axis data to parse
    lowest : float
        Lowest value to plot on this axis
    highest : float
        Highest value to plot on this axis

    RETURNS
    -------
    axis_data : an AxisData object
        The AxisData object that contains the resulting axis data
    """

    if axis_data is None:
        return EmptyAxisData()

    else:
        try:
            column = int(axis_data[0])

            if len(axis_data) > 2:
                raise ValueError
        except ValueError:
            print("Warning: axis data not recognized: {}".format(axis_data))
            return EmptyAxisData()

        if len(axis_data) == 1:
            return UnlabelledAxisData(column, lowest=lowest,
                                      highest=highest)

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

            return AxisData(label=label, column=column,
                            lowest=lowest, highest=highest)


def interactive_plot(pele_data, plot_data):
    """
    It launches an interactive plot.

    Parameters
    ----------
    pele_data : a Pandas.DataFrame object
        The DataFrame that contains the PELE simulation data
    plot_data : a PlotData object
        The PlotData containing the information to plot in each
        axes
    """
    from matplotlib import pyplot

    x_values = plot_data.get_xs_from_pele_data(pele_data)
    y_values = plot_data.get_ys_from_pele_data(pele_data)
    z_values = plot_data.get_zs_from_pele_data(pele_data)

    if z_values is not None:
        if plot_data.z_data.highest is None:
            z_max = max(z_values)
        else:
            z_max = plot_data.z_data.highest

        if plot_data.z_data.lowest is None:
            z_min = min(z_values)
        else:
            z_min = plot_data.z_data.lowest

        cmap = pyplot.cm.plasma
        norm = pyplot.Normalize(z_min, z_max)
    else:
        cmap = None
        norm = None

    fig, ax = pyplot.subplots()

    scatter = pyplot.scatter(x_values, y_values, c=z_values, cmap=cmap,
                             norm=norm)

    ax.margins(0.05)
    ax.set_facecolor('lightgray')
    pyplot.ylabel(plot_data.y_data.label)
    pyplot.xlabel(plot_data.x_data.label)
    pyplot.xlim([plot_data.x_data.lowest, plot_data.x_data.highest])
    pyplot.ylim([plot_data.y_data.lowest, plot_data.y_data.highest])

    annotation = ax.annotate("", xy=(0, 0), xytext=(20, 20),
                             textcoords="offset points",
                             bbox=dict(boxstyle="round", fc="w"),
                             arrowprops=dict(arrowstyle="->"))
    annotation.set_visible(False)

    # Activate the colorbar only if the Z axis contains data to plot
    if not isinstance(plot_data.z_data, EmptyAxisData):
        cbar = pyplot.colorbar(scatter, drawedges=False)
        cbar.ax.set_ylabel(plot_data.z_data.label)

    def update_annot(idx):
        """Update the information box of the selected point"""
        position = scatter.get_offsets()[idx["ind"][0]]
        annotation.xy = position

        # Extract point information
        pele_data_row = pele_data.iloc[idx["ind"][0], :]
        epoch = str(pele_data_row['epoch'])
        trajectory = os.path.basename(pele_data_row['trajectory'])
        model = str(pele_data_row['numberOfAcceptedPeleSteps'] + 1)

        # Annotate information
        annotation.set_text(("Epoch: " + epoch + "\n" +
                             "Trajectory: " + trajectory + "\n" +
                             "Model: " + model))

        # Set annotation's face color
        if cmap is not None:
            color = cmap(norm(z_values[idx["ind"][0]]))
            annotation.get_bbox_patch().set_facecolor(color)

    def hover(event):
        """Action to perform when hovering the mouse on a point"""
        visible = annotation.get_visible()
        if event.inaxes == ax:
            inside, idx = scatter.contains(event)
            if inside:
                update_annot(idx)
                annotation.set_visible(True)
                fig.canvas.draw_idle()
            else:
                if visible:
                    annotation.set_visible(False)
                    fig.canvas.draw_idle()

    # Respond to mouse motion
    fig.canvas.mpl_connect("motion_notify_event", hover)

    pyplot.show()


def scatter_plot(pele_data, plot_data):
    """
    It launches an scatter plot.

    Parameters
    ----------
    pele_data : a Pandas.DataFrame object
        The DataFrame that contains the PELE simulation data
    plot_data : a PlotData object
        The PlotData containing the information to plot in each
        axes
    """
    from matplotlib import pyplot

    x_values = plot_data.get_xs_from_pele_data(pele_data)
    y_values = plot_data.get_ys_from_pele_data(pele_data)
    z_values = plot_data.get_zs_from_pele_data(pele_data)

    fig, ax = pyplot.subplots()

    if z_values is None:
        ax.scatter(x_values, y_values, s=20)

    else:
        scatter = ax.scatter(x_values, y_values, c=z_values, s=20)
        cbar = pyplot.colorbar(scatter)
        cbar.ax.set_ylabel(plot_data.z_data.label)

    ax.margins(0.05)
    ax.set_facecolor('lightgray')
    pyplot.ylabel(plot_data.y_data.label)
    pyplot.xlabel(plot_data.x_data.label)
    pyplot.xlim([plot_data.x_data.lowest, plot_data.x_data.highest])
    pyplot.ylim([plot_data.y_data.lowest, plot_data.y_data.highest])

    pyplot.show()


# Main workflow to be executed
if __name__ == "__main__":
    # Import external libraries
    import pandas as pd

    # Avoid Python2
    from .Checker.python_version import check_python_version

    check_python_version()

    # Parse command-line arguments
    csv_file, results_folder, output_folder, report_name, \
        trajectory_name, plot_type, xdata, ydata, zdata,\
        xlowest, xhighest, ylowest, yhighest, zlowest, zhighest, = \
        parse_args()

    # Print header
    from pele_platform.constants import constants
    print(constants.plotter_version_header.format(__version__))

    print_parameters(csv_file, results_folder, output_folder,
                     report_name, trajectory_name, plot_type,
                     xdata, ydata, zdata)

    if csv_file is not None:
        pele_data = pd.read_csv(csv_file)

    elif results_folder is not None:
        pele_data = pd.read_csv(os.path.join(results_folder, 'data.csv'))

    elif output_folder is not None:
        from pele_platform.analysis import DataHandler

        data_handler = DataHandler(sim_path=output_folder,
                                   report_name=report_name,
                                   trajectory_name=trajectory_name,
                                   skip_initial_structures=False)
        pele_data = data_handler.get_reports_dataframe()

    else:
        raise Exception('Data is missing, check your arguments.')  # This message should never be prompted

    x_data = parse_axis_data(xdata, xlowest, xhighest)
    y_data = parse_axis_data(ydata, ylowest, yhighest)
    z_data = parse_axis_data(zdata, zlowest, zhighest)

    plot_data = PlotData(x_data, y_data, z_data)

    if plot_type is None:
        plot_type = DEFAULT_PLOT_TYPE

    if plot_type.lower() == 'interactive':
        interactive_plot(pele_data, plot_data)
    elif plot_type.lower() == 'scatter':
        scatter_plot(pele_data, plot_data)




