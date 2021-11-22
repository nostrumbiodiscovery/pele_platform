# -*- coding: utf-8 -*-
"""
This is the plotter module of the platform that can take regular PELE
outputs and display them on a variety of graphs.
"""

__author__ = "Nostrum Biodiscovery"
__email__ = "pelesupport@nostrumbiodiscovery.com"
__license__ = "Apache-2.0"
__version__ = "1.1"


# General imports
import os


# Constants
PLOT_TYPE_CHOICES = ('scatter', 'interactive', 'density', 'histogram')
DEFAULT_PLOT_TYPE = 'scatter'
COLORS = {'blue': ('lightskyblue', 'royalblue'),
          'red': ('#f5bcbc', 'firebrick'),
          'green': ('palegreen', 'seagreen'),
          'purple': ('#d9c6ec', '#8000ff'),
          'orange': ('#facc9e', '#ff6600')}


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

    def get_xs_from_pele_data(self, pele_data, sanitize=True):
        """
        It returns the X values to plot from a DataFrame containing
        PELE data.

        Parameters
        ----------
        pele_data : a Pandas.DataFrame object
            The DataFrame containing the PELE simulation data from
            where the values will be extracted
        sanitize : bool
            Whether to sanitize data or not. Default is 'True'

        Returns
        -------
        xs : a Pandas.DataFrame object
            The DataFrame containing the values to plot on the X axis
        """
        if self.x_data.category == 'report_column':
            xs = pele_data.iloc[:, self.x_data.column - 1]

            if sanitize:
                import pandas as pd
                pd.to_numeric(xs, errors='coerce').fillna(-1)

            return xs
        else:
            raise NotImplementedError()

    def get_ys_from_pele_data(self, pele_data, sanitize=True):
        """
        It returns the Y values to plot from a DataFrame containing
        PELE data.

        Parameters
        ----------
        pele_data : a Pandas.DataFrame object
            The DataFrame containing the PELE simulation data from
            where the values will be extracted
        sanitize : bool
            Whether to sanitize data or not. Default is 'True'

        Returns
        -------
        ys : a Pandas.DataFrame object
            The DataFrame containing the values to plot on the Y axis
        """
        if self.y_data.category == 'report_column':
            ys = pele_data.iloc[:, self.y_data.column - 1]

            if sanitize:
                import pandas as pd
                pd.to_numeric(ys, errors='coerce').fillna(-1)

            return ys
        else:
            raise NotImplementedError()

    def get_zs_from_pele_data(self, pele_data, sanitize=True):
        """
        It returns the Z values to plot from a DataFrame containing
        PELE data.

        Parameters
        ----------
        pele_data : a Pandas.DataFrame object
            The DataFrame containing the PELE simulation data from
            where the values will be extracted
        sanitize : bool
            Whether to sanitize data or not. Default is 'True'

        Returns
        -------
        zs : a Pandas.DataFrame object
            The DataFrame containing the values to plot on the Z axis
        """
        if isinstance(self.z_data, EmptyAxisData):
            return None

        if self.z_data.category == 'report_column':
            zs = pele_data.iloc[:, self.z_data.column - 1]

            if sanitize:
                import pandas as pd
                zs = pd.to_numeric(zs, errors='coerce').fillna(-1)

            return zs
        else:
            raise NotImplementedError()

    def is_plottable(self):
        """
        It notifies if the current plot data contains enough data to
        generate a plot.

        Returns
        -------
        answer : bool
            Whether this PlotData object can be visualized in a plot
            or not
        """
        answer = (not isinstance(self.x_data, EmptyAxisData) and
                  not isinstance(self.y_data, EmptyAxisData))

        return answer


class PlotAppearance(object):
    """
    It handles the appearance settings of the plot.
    """

    def __init__(self, colormap_name, plot_color, background_color,
                 display_edges=None, n_levels=5, n_bins=[100, 100],
                 lines=[], title=None, hide_logo=False):
        """
        It initializes a PlotAppearance object.

        Parameters
        ----------
        colormap_name : str
            The name of the colormap to represent the values
            of the Z axis
        plot_color : str
            The color for the plot. One of ['blue', 'red', 'green',
            'purple', 'orange']
        background_color : str
            The background color for the plot. Default is 'None'
        display_edges : bool
            Whether to display level edges or not in a density plot.
            Default is False
        n_levels : int
            Number of levels to display in a density plot. Default is 5
        n_bins : tuple[int, int]
            Number of bins to display in the histogram plot, first element
            corresponds to the X axis and the second to the Y axis
        lines : list[Line]
            A list of Line objects. Default is an empty list
        title : str
            The title to print in the plot. Default is None
        hide_logo : bool
            Whether to hide NBD logo on the plot or not. Default is False
        """
        self.colormap_name = colormap_name
        self.plot_color = plot_color
        self.background_color = background_color
        self.display_edges = display_edges
        self.n_levels = n_levels
        self.n_bins = n_bins
        self.lines = lines
        self.title = title
        self.hide_logo = hide_logo


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


class Line(object):
    """
    It represents a line to be added in the plot.
    """

    def __init__(self, intercept, color=None, vertical=True):
        """
        It initializes a Line object.

        Parameters
        ----------
        intercept : float
           The value where the line intercepts the axis
        color : str
            The color name for this line. It must be a CSS4 compatible
            color name. Default is None, which will set color red
        vertical : bool
           Whether it is a vertical or an horizontal line. Default is
           'True'
        """

        self.intercept = intercept
        if color is None:
            self.color = 'red'
        else:
            self.color = color
        self.vertical = vertical

    def __str__(self):
        """
        It returns the string representation of this Line instance.

        Returns
        -------
        str_representation : str
            The string representation of this Line instance
        """
        str_representation = f"{self.intercept} {self.color}"

        return str_representation

    def __repr__(self):
        """
        It returns the representation of this Line instance.

        Returns
        -------
        repr : str
            The representation of this Line instance
        """
        repr = str(self)

        return repr


class Filter(object):
    """
    It represents a filter to be applied to PELE data.
    """

    def __init__(self, column, condition, threshold):
        """
        It initializes a Filter object.

        Parameters
        ----------
        column : int
           The column number corresponding to the metric where
           the filter will be applied
        condition : str
            It represents the condition to apply in the filtering.
            One of ['<', '==', '>', '<=', '>=']
        threshold : float
           Threshold value to be applied to the column data.
        """

        if not isinstance(column, int):
            raise TypeError('Wrong type for column. Integer is expected, ' +
                            f'got {column}')

        if condition not in ['<', '==', '>', '<=', '>=']:
            raise NameError('Wrong condition. One of ' +
                            '[\'<\', \'=\', \'>\', \'<=\', \'>=\'] ' +
                            f'is expected, got {condition}')

        if not isinstance(threshold, (float, int)):
            raise TypeError('Wrong type for column. Float is expected, ' +
                            f'got {threshold}')

        self.column = column
        self.condition = condition
        self.threshold = threshold

    def __str__(self):
        """
        It sets the string representation of the Filter class.

        Returns
        -------
        string_repr : str
            The string representation of a Filter object
        """

        return f'{self.column} {self.condition} {self.threshold}'

    def apply(self, pele_data):
        """
        It applies the filter to the PELE data that is supplied.

        Parameters
        ----------
        pele_data : a Pandas.DataFrame object
            The DataFrame that contains the PELE simulation data

        Returns
        -------
        filtered_pele_data : a Pandas.DataFrame object
            The DataFrame that contains the filtered data
        """
        initial_points = len(pele_data)

        try:
            pele_data.iloc[:, self.column - 1] = \
                pele_data.iloc[:, self.column - 1] .astype(float)
        except ValueError:
            print(f' - Values in column {self.column} are not numeric ' +
                  'and cannot be filtered')
            return pele_data

        if self.condition == '<':
            pele_data = pele_data[pele_data.iloc[:, self.column - 1]
                                  < self.threshold]
        elif self.condition == '>':
            pele_data = pele_data[pele_data.iloc[:, self.column - 1]
                                  > self.threshold]
        elif self.condition == '==':
            pele_data = pele_data[pele_data.iloc[:, self.column - 1]
                                  == self.threshold]
        elif self.condition == '>=':
            pele_data = pele_data[pele_data.iloc[:, self.column - 1]
                                  >= self.threshold]
        elif self.condition == '<=':
            pele_data = pele_data[pele_data.iloc[:, self.column - 1]
                                  <= self.threshold]

        final_points = len(pele_data)

        print(f' - {initial_points} points were reduced to {final_points}')

        return pele_data


# Methods
def parse_line_data(line_data, vertical):
    """
    It parses the supplied line data and initializes the corresponding
    Line objects.

    Parameters
    ----------
    line_data : list[tuple[str, str]] or list[tuple[str]]
        It is a list of 1 or 2 dimensional tuples containing the
        intersection of each line and, if it is also supplied, its
        color. First tuple element is the intersection, the second
        element is not mandatory and specifies the custom color
        for that line
    vertical : bool
        Whether line data belongs to vertical or horizontal lines

    Returns
    -------
    lines : list[Line]
        It is a list of the resulting Line objects
    """
    class LineParserException(BaseException):
        """It sets a line parser exception."""
        pass

    lines = list()

    for one_line_data in line_data:
        try:
            if not isinstance(one_line_data, (list, tuple)):
                raise LineParserException()

            intercept = float(one_line_data[0])

            if len(one_line_data) == 1:
                color = None

            elif len(one_line_data) == 2:
                color = str(one_line_data[1])

                from matplotlib import colors as mcolors
                if color not in mcolors.CSS4_COLORS:
                    raise ValueError('wrong color: ' +
                                     f'{color}, it must be a ' +
                                     'CSS4-compatible color name.')

            else:
                raise LineParserException('wrong format, line input can ' +
                                          'only contain either 1 or 2 fields')

            line = Line(intercept=intercept, color=color,
                        vertical=vertical)
            lines.append(line)

        except (LineParserException, ValueError) as error:
            print(f"Warning: line data not recognized: {one_line_data},",
                  str(error))

    return lines


def parse_filters(filters_data):
    """
    It parses the supplied filter data and initializes the corresponding
    Filter objects.

    Parameters
    ----------
    filters_data : list[tuple[int, str, float]]
        It is a list 3 dimensional tuples containing the filters
        to be applied into the PELE data

    Returns
    -------
    filters : list[Filter]
        It is a list of parsed filters
    """
    class FilterParserException(BaseException):
        """It sets a line parser exception."""
        pass

    filters = list()

    for one_filter_data in filters_data:
        try:
            if not isinstance(one_filter_data, (list, tuple)):
                raise FilterParserException()

            if len(one_filter_data) != 3:
                raise FilterParserException('Invalid filter format. ' +
                                            'Expected: \'column\', '+
                                            '\'condition\', ' +
                                            '\'threshold\'], ' +
                                            f'got {one_filter_data}')

            try:
                column = int(one_filter_data[0])
            except TypeError:
                raise FilterParserException('Invalid column. ' +
                                            'Expected an integer, got ' +
                                            f'{one_filter_data[0]}')

            condition = one_filter_data[1]

            if condition.lower() == 'eq':
                condition = '=='
            elif condition.lower() == 'lt':
                condition = '<'
            elif condition.lower() == 'gt':
                condition = '>'
            elif condition.lower() == 'le':
                condition = '<='
            elif condition.lower() == 'ge':
                condition = '>='

            if condition not in ['<', '==', '>', '<=', '>=']:
                raise FilterParserException('Wrong condition. One of ' +
                                            '[\'<\', \'==\', \'>\', ' +
                                            '\'<=\', \'>=\', \'eq\', ' +
                                            '\'lt\', \'eq\', \'gt\', ' +
                                            '\'le\', \'ge\'] '+
                                            'is expected, ' +
                                            f'got {condition}')

            try:
                threshold = float(one_filter_data[2])
            except TypeError:
                raise FilterParserException('Invalid column. ' +
                                            'Expected a float, got ' +
                                            f'{one_filter_data[0]}')

            filter = Filter(column=column, condition=condition,
                            threshold=threshold)
            filters.append(filter)

        except FilterParserException as error:
            print(f"Warning: filter data not recognized: {one_filter_data}, ",
                  str(error))

    return filters


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
    xlowest : float
        Lowest value to plot on the X axis
    xhighest : float
        Highest value to plot on the X axis
    ylowest : float
        Lowest value to plot on the Y axis
    yhighest : float
        Highest value to plot on the Y axis
    zlowest : float
        Lowest value to plot on the Z axis
    zhighest : float
        Highest value to plot on the Z axis
    colormap : str
        The name of the colormap to represent the values
        of the Z axis
    color : str
        The color for the plot. One of ['blue', 'red', 'green',
        'purple', 'orange']
    with_edges : bool
        Display edges of levels in the density plot
    n_levels : int
        Number of levels to display in the density plot when
        edges are shown
    n_bins : tuple[int, int]
        Number of bins to display in the histogram plot, first element
        corresponds to the X axis and the second to the Y axis
    background_color : str
        The background color for the plot. It must be a CSS4 compatible
        color name
    vertical_lines : list[Line]
        A list of all vertical Line objects
    horizontal_lines : list[Line]
        A list of all horizontal Line objects
    filters : list[Filter]
        The list of filters to apply to PELE data
    title : str
        The title to print in the plot
    hide_logo : bool
        Whether to hide NBD logo on the plot or not
    save_to : str
        If it is not None, the resulting plot will be saved to this path
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
                        choices=PLOT_TYPE_CHOICES, default=DEFAULT_PLOT_TYPE,
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
    parser.add_argument('--colormap', type=str, default='gnuplot',
                        help='The colormap to employ to represent ' +
                        'the values of the Z axis')
    parser.add_argument('--color', type=str,
                        choices=COLORS, default='blue',
                        help='Plot color. One of ' + str(COLORS))
    parser.add_argument('--with_edges', dest='with_edges',
                        action='store_true',
                        help='Display edges of levels in the density plot')
    parser.add_argument('--n_levels', type=int, default=5,
                        help='Number of levels to display in the density ' +
                        'plot when edges are shown')
    parser.add_argument('--n_bins', type=int, default=[100, ], nargs='*',
                        help='Number of bins to use in a histogram plot. ' +
                        'Either one value to apply to both axis or two ' +
                        'values corresponding to X and Y axes, respectively.')
    parser.add_argument('--background_color', type=str, default='white',
                        help='Background color for the plot.')
    parser.add_argument('--vertical_line', action='append', type=str,
                        help='It adds a vertical line to the intercept ' +
                        'that is supplied.', nargs='*',
                        metavar=('INT', 'STR'))
    parser.add_argument('--horizontal_line', action='append', type=str,
                        help='It adds an horizontal line to the intercept ' +
                        'that is supplied.', nargs='*',
                        metavar=('INT', 'STR'))
    parser.add_argument('--filter', action='append', type=str,
                        help='It applies a filter to the PELE dataset. ' +
                        'It expects the following data: metric_column + ' +
                        'condition + threshold_value. Condition must be ' +
                        'one of [\'<\', \'>\', \'==\', \'<=\', \'>=\', ' +
                        '\'lt\', \'eq\', \'gt\', \'le\', \'ge\'].',
                        nargs=3, metavar=('INT', 'STR', 'FLOAT'))
    parser.add_argument('--title', type=str, default=None,
                        help='The title to print in the plot.')
    parser.add_argument('--hide_logo', dest='hide_logo',
                        action='store_true',
                        help='It hides NBD logo')
    parser.add_argument('-s', '--save_to', type=str,
                       default=None, help='Instead of displaying the ' +
                       'plot, it will save it to the supplied path.')
    # Set defaults
    parser.set_defaults(with_edges=False)
    parser.set_defaults(hide_logo=False)

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
    colormap = parsed_args.colormap
    color = parsed_args.color
    with_edges = parsed_args.with_edges
    n_levels = parsed_args.n_levels
    n_bins = parsed_args.n_bins
    background_color = parsed_args.background_color
    vertical_lines = parsed_args.vertical_line
    horizontal_lines = parsed_args.horizontal_line
    filters = parsed_args.filter
    title = parsed_args.title
    hide_logo = parsed_args.hide_logo
    save_to = parsed_args.save_to

    if vertical_lines is None:
        vertical_lines = list()
    else:
        vertical_lines = parse_line_data(vertical_lines, vertical=True)

    if horizontal_lines is None:
        horizontal_lines = list()
    else:
        horizontal_lines = parse_line_data(horizontal_lines, vertical=False)

    if filters is None:
        filters = list()
    else:
        filters = parse_filters(filters)

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

    from matplotlib import colors as mcolors
    if (background_color is not None and
            background_color not in mcolors.CSS4_COLORS):
        raise ValueError(f'Wrong background color: {background_color}, ' +
                         'it must be a CSS4-compatible color name.')

    if save_to is not None:
        if os.path.isdir(save_to):
            plot_name = f"pele_{plot_type.lower()}_plot.png"
            save_to = os.path.join(save_to, plot_name)
        else:
            if (os.path.dirname(save_to) != '' and
                    not os.path.isdir(os.path.dirname(save_to))):
                raise ValueError(f'Plot cannot be saved at {save_to}. '
                                 'Directory does not exist.')
        if os.path.splitext(os.path.basename(save_to))[1] == '':
            save_to = save_to + '.png'
        if os.path.splitext(os.path.basename(save_to))[1] != '.png':
            raise ValueError('Plot can only be saved with PNG format')

    if len(n_bins) == 1:
        n_bins = n_bins * 2
    elif len(n_bins) != 2:
        raise ValueError('Invalid number of bins: either one value to ' +
                         'apply it to both axis or two values corresponding ' +
                         'to X and Y axes, respectively, are accepted')

    return csv_file, results_folder, output_folder, \
        report_name, trajectory_name, plot_type, xdata, ydata, zdata, \
        xlowest, xhighest, ylowest, yhighest, zlowest, zhighest, \
        colormap, color, with_edges, n_levels, n_bins, \
        background_color, vertical_lines, horizontal_lines, filters, \
        title, hide_logo, save_to


def print_parameters(csv_file, results_folder, output_folder,
                     report_name, trajectory_name,
                     plot_type, xdata, ydata, zdata,
                     xlowest, xhighest, ylowest, yhighest,
                     zlowest, zhighest, colormap, color,
                     with_edges, n_levels, n_bins, background_color,
                     vertical_lines, horizontal_lines, filters, title,
                     hide_logo, save_to):
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
    xlowest : float
        Lowest value to plot on the X axis
    xhighest : float
        Highest value to plot on the X axis
    ylowest : float
        Lowest value to plot on the Y axis
    yhighest : float
        Highest value to plot on the Y axis
    zlowest : float
        Lowest value to plot on the Z axis
    zhighest : float
        Highest value to plot on the Z axis
    colormap : str
        The name of the colormap to represent the values
        of the Z axis
    color : str
        The color for the plot. One of ['blue', 'red', 'green',
        'purple', 'orange']
    with_edges : bool
        Display edges of levels in the density plot
    n_levels : int
        Number of levels to display in the density plot when
        edges are shown
    n_bins : tuple[int, int]
        Number of bins to display in the histogram plot, first element
        corresponds to the X axis and the second to the Y axis
    background_color : str
        The background color for the plot. It must be a CSS4 compatible
        color name
    vertical_lines : list[Line]
        A list of all vertical Line objects
    horizontal_lines : list[Line]
        A list of all horizontal Line objects
    filters : list[Filter]
        The list of filters to apply to PELE data
    title : str
        The title to print in the plot
    hide_logo : bool
        Whether to hide NBD logo on the plot or not
    save_to : str
        If it is not None, the resulting plot will be saved to this path
    """

    max_len = max([len(str(item)) if item is not None else 1
                   for item in [csv_file, results_folder,
                                output_folder, report_name,
                                trajectory_name, plot_type,
                                xdata, ydata, zdata,
                                xlowest, xhighest,
                                ylowest, yhighest,
                                zlowest, zhighest,
                                colormap, color, with_edges,
                                n_levels, n_bins, background_color,
                                vertical_lines, horizontal_lines,
                                title, hide_logo, save_to]])

    print('-> Input parameters:')
    print(' - csv_file:         ',
          ' ' * (max_len - (len(str(csv_file))
                            if csv_file is not None else 1)),
          '-' if csv_file is None else csv_file)
    print(' - results_folder:   ',
          ' ' * (max_len - (len(str(results_folder))
                            if results_folder is not None else 1)),
          '-' if results_folder is None else results_folder)
    print(' - output_folder:    ',
          ' ' * (max_len - (len(str(output_folder))
                            if output_folder is not None else 1)),
          '-' if output_folder is None else output_folder)
    print(' - report_name:      ',
          ' ' * (max_len - (len(str(report_name))
                            if report_name is not None else 1)),
          '-' if report_name is None else report_name)
    print(' - trajectory_name:  ',
          ' ' * (max_len - (len(str(trajectory_name))
                            if trajectory_name is not None else 1)),
          '-' if trajectory_name is None else trajectory_name)
    print(' - plot_type:        ',
          ' ' * (max_len - (len(str(plot_type))
                            if plot_type is not None else 1)),
          '-' if plot_type is None else plot_type)
    print(' - xdata:            ',
          ' ' * (max_len - (len(str(xdata))
                            if xdata is not None else 1)),
          '-' if xdata is None else xdata)
    print(' - ydata:            ',
          ' ' * (max_len - (len(str(ydata))
                            if ydata is not None else 1)),
          '-' if ydata is None else ydata)
    print(' - zdata:            ',
          ' ' * (max_len - (len(str(zdata))
                            if zdata is not None else 1)),
          '-' if zdata is None else zdata)
    print(' - xlowest:          ',
          ' ' * (max_len - (len(str(xlowest))
                            if xlowest is not None else 1)),
          '-' if xlowest is None else xlowest)
    print(' - xhighest:         ',
          ' ' * (max_len - (len(str(xhighest))
                            if xhighest is not None else 1)),
          '-' if xhighest is None else xhighest)
    print(' - ylowest:          ',
          ' ' * (max_len - (len(str(ylowest))
                            if ylowest is not None else 1)),
          '-' if ylowest is None else ylowest)
    print(' - yhighest:         ',
          ' ' * (max_len - (len(str(yhighest))
                            if yhighest is not None else 1)),
          '-' if yhighest is None else yhighest)
    print(' - zlowest:          ',
          ' ' * (max_len - (len(str(zlowest))
                            if zlowest is not None else 1)),
          '-' if zlowest is None else zlowest)
    print(' - zhighest:         ',
          ' ' * (max_len - (len(str(zhighest))
                            if zhighest is not None else 1)),
          '-' if zhighest is None else zhighest)
    print(' - colormap:         ',
          ' ' * (max_len - (len(str(colormap))
                            if colormap is not None else 1)),
          '-' if colormap is None else colormap)
    print(' - color:            ',
          ' ' * (max_len - (len(str(color))
                            if color is not None else 1)),
          '-' if color is None else color)
    print(' - with_edges:       ',
          ' ' * (max_len - (len(str(with_edges))
                            if with_edges is not None else 1)),
          '-' if with_edges is None else with_edges)
    print(' - n_levels:         ',
          ' ' * (max_len - (len(str(n_levels))
                            if n_levels is not None else 1)),
          '-' if n_levels is None else n_levels)
    print(' - n_bins:           ',
          ' ' * (max_len - (len(str(n_bins))
                            if n_bins is not None else 1)),
          '-' if n_bins is None else n_bins)
    print(' - background_color: ',
          ' ' * (max_len - (len(str(background_color))
                            if background_color is not None else 1)),
          '-' if background_color is None else background_color)
    print(' - vertical_lines:   ',
          ' ' * (max_len - (len(str(vertical_lines))
                            if vertical_lines is not None else 1)),
          '-' if vertical_lines is None else vertical_lines)
    print(' - horizontal_lines: ',
          ' ' * (max_len - (len(str(horizontal_lines))
                            if horizontal_lines is not None else 1)),
          '-' if horizontal_lines is None else horizontal_lines)
    if len(filters) == 0:
        print(' - filters:          ',
              ' ' * (max_len - 1),
              '-')
    else:
        for i, filter in enumerate(filters, start=1):
            print(f' - filter{i}: ',
                  ' ' * (9 - len(str(i))),
                  ' ' * (max_len - (len(str(filter)))),
                  str(filter))
    print(' - title:            ',
          ' ' * (max_len - (len(str(title))
                            if title is not None else 1)),
          '-' if title is None else title)
    print(' - hide_logo:        ',
          ' ' * (max_len - (len(str(hide_logo))
                            if hide_logo is not None else 1)),
          '-' if hide_logo is None else hide_logo)
    print(' - save_to:          ',
          ' ' * (max_len - (len(str(save_to))
                            if save_to is not None else 1)),
          '-' if save_to is None else save_to)
    print()


def request_axis_data(axis_name, default_answer, pele_data, optional=False):
    """
    It asks the user for the axis data.

    Parameters
    ----------
    axis_name : str
        The name of the axis whose data is requested
    default_answer : str
        The default answer to use. One of ['y', 'n']
    pele_data : a Pandas.DataFrame object
        The DataFrame that contains the PELE simulation data
    optional : bool
        Whether the user will be forced to supply axis data or they
        will be able to skip the request. Default is 'False'

    Returns
    -------
    axis_data : an AxisData object
        The AxisData object that contains the resulting axis data
    """
    assert default_answer.lower() in ['y', 'n'], \
        'Wrong default answer supplied'

    if optional:
        optional_str = ' (optional)'

        if default_answer == 'y':
            default_opt = '[y]/n'
        else:
            default_opt = 'y/[n]'

    else:
        optional_str = ''
        default_opt = ''

    print(f'-> No data to plot on axis {axis_name} was supplied' +
          f'{optional_str}.')

    if optional:
        answer = input(f'  Would you like to supply it now ' +
                       f'({default_opt})? ')

        if answer == '':
            answer = default_answer

        if answer.lower() == 'n':
            print()
            return EmptyAxisData()

        if answer.lower() not in ['y', 'n']:
            print(f' Unexpected answer.')
            print()
            return request_axis_data(axis_name, default_answer,
                                     pele_data)

    report_columns = list(pele_data.columns)
    max_len = max([len(metric) for metric in report_columns])
    char_counter = 2
    final_break = False
    print(f' - Available metrics:')
    for idx, metric in enumerate(report_columns, start=1):
        if not final_break:
            print('  ', end='')
        to_print = '{:2d}. {}{}'.format(idx,
                                        metric,
                                        ' ' * (max_len - len(metric) + 1))
        if char_counter + len(to_print) > 80:
            print(to_print)
            char_counter = 2
            final_break = False
        else:
            print(to_print, end='')
            char_counter += len(to_print)
            final_break = True

    if final_break:
        print()

    column = input(f' - Metric column index: ')

    if not column.isdigit() or not float(column).is_integer():
        print(f' Unexpected answer: column index must be an integer.')
        print()
        return request_axis_data(axis_name, default_answer,
                                 pele_data)

    column = int(column)

    if column < 1 or column > len(report_columns):
        print(f' Unexpected answer: column index must be inside ' +
              'range: [1, {}].'.format(len(report_columns)))
        print()
        return request_axis_data(axis_name, default_answer,
                                 pele_data)

    label = input(f' - Axis label [optional]: ')
    if label == '':
        label = report_columns[column - 1]

    label = add_units(label)

    print()
    return AxisData(label, column)


def add_units(metric_name):
    """
    Given a metric name, it adds the units according to its input.

    Parameters
    ----------
    metric_name : str
        The name of the metric where the units will be added

    Returns
    -------
    label : str
        The new label of the metric containing the units that
        were inferred from its original content
    """

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

    return label


def get_colormap(colormap_name):
    """
    Given a colormap name, it returns the corresponding colormap object.

    Parameters
    ----------
    colormap_name : str
        The name of the colormap

    Returns
    -------
    cmap : a Matplotlib.cm.Colormap object
        The requested colormap object
    """
    from matplotlib import pyplot
    from matplotlib.colors import ListedColormap
    import numpy as np

    if colormap_name.lower() == 'plasma':
        cmap = pyplot.cm.plasma
    elif colormap_name.lower() == 'magma':
        cmap = pyplot.cm.magma
    elif colormap_name.lower() == 'turbo':
        cmap = pyplot.cm.turbo
    elif colormap_name.lower() == 'jet':
        cmap = pyplot.cm.jet
    elif colormap_name.lower() == 'gnuplot':
        cmap = pyplot.cm.gnuplot
    elif colormap_name.lower() == 'gnuplot2':
        cmap = pyplot.cm.gnuplot2
    elif colormap_name.lower() == 'nipy_spectral':
        cmap = pyplot.cm.nipy_spectral
    elif colormap_name.lower() == 'spectral':
        cmap = pyplot.cm.Spectral
    elif colormap_name.lower() == 'cividis':
        cmap = pyplot.cm.cividis
    elif colormap_name.lower() == 'inferno':
        cmap = pyplot.cm.inferno
    elif colormap_name.lower() == 'autumn':
        cmap = pyplot.cm.autumn
    elif colormap_name.lower() == 'winter':
        cmap = pyplot.cm.winter
    elif colormap_name.lower() == 'spring':
        cmap = pyplot.cm.spring
    elif colormap_name.lower() == 'summer':
        cmap = pyplot.cm.summer
    elif colormap_name.lower() == 'wistia':
        cmap = pyplot.cm.Wistia
    elif colormap_name.lower() == 'copper':
        cmap = pyplot.cm.copper
    elif colormap_name.lower() == 'blues':
        cmap = pyplot.cm.Blues
    elif colormap_name.lower() == 'reducedblues':
        blues = pyplot.cm.get_cmap('Blues', 512)
        cmap = ListedColormap(blues(np.linspace(0.95, 0.4, 256)))
    elif colormap_name.lower() == 'reducedgreens':
        greens = pyplot.cm.get_cmap('Greens', 512)
        cmap = ListedColormap(greens(np.linspace(0.95, 0.4, 256)))
    elif colormap_name.lower() == 'reducedreds':
        reds = pyplot.cm.get_cmap('Reds', 512)
        cmap = ListedColormap(reds(np.linspace(0.95, 0.4, 256)))
    elif colormap_name.lower() == 'reducedpurples':
        purples = pyplot.cm.get_cmap('Purples', 512)
        cmap = ListedColormap(purples(np.linspace(0.95, 0.4, 256)))
    elif colormap_name.lower() == 'reducedoranges':
        oranges = pyplot.cm.get_cmap('Oranges', 512)
        cmap = ListedColormap(oranges(np.linspace(0.95, 0.4, 256)))
    else:
        try:
            cmap = pyplot.cm.get_cmap(colormap_name)
        except ValueError:
            raise NameError(f'Unknown colormap name: \'{colormap_name}\'')

    return cmap


def parse_axis_data(axis_data, lowest, highest, pele_data):
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
    pele_data : a Pandas.DataFrame object
        The DataFrame that contains the PELE simulation data

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
            report_columns = list(pele_data.columns)
            metric_name = report_columns[column - 1]

        else:
            metric_name = axis_data[1]

        label = add_units(metric_name)

        return AxisData(label=label, column=column,
                        lowest=lowest, highest=highest)


def interactive_plot(pele_data, plot_data, plot_appearance, save_to):
    """
    It generates an interactive plot.

    Parameters
    ----------
    pele_data : a Pandas.DataFrame object
        The DataFrame that contains the PELE simulation data
    plot_data : a PlotData object
        The PlotData containing the information to plot in each
        axes
    plot_appearance : a PlotAppearance object
        The appearance settings for the plot
    save_to : str
        If it is not None, the resulting plot will be saved to this path
    """
    from matplotlib import pyplot

    # Check parameters
    if save_to is not None:
        raise ValueError('Interactive plots cannot be saved.')

    # Set seaborn style
    pyplot.style.use("seaborn-ticks")

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

        cmap = get_colormap(plot_appearance.colormap_name)
        norm = pyplot.Normalize(z_min, z_max)
        color = None
    else:
        cmap = None
        norm = None
        color = COLORS[plot_appearance.plot_color][1]

    fig, ax = pyplot.subplots()

    scatter = pyplot.scatter(x_values, y_values, c=z_values, cmap=cmap,
                             norm=norm, color=color)

    # Axes settings
    title = plot_appearance.title
    if title is None:
        title = "PELE Interactive Plot"
    ax.set_title(title, fontweight='bold', pad=20)
    ax.margins(0.05)
    ax.set_facecolor(plot_appearance.background_color)
    pyplot.ylabel(plot_data.y_data.label, fontweight='bold')
    pyplot.xlabel(plot_data.x_data.label, fontweight='bold')
    pyplot.xlim([plot_data.x_data.lowest, plot_data.x_data.highest])
    pyplot.ylim([plot_data.y_data.lowest, plot_data.y_data.highest])

    annotation = ax.annotate("", xy=(0, 0), xytext=(20, 20),
                             textcoords="offset points",
                             bbox=dict(boxstyle="round", fc="w"),
                             arrowprops=dict(arrowstyle="->"))
    annotation.set_visible(False)

    # Add lines
    for line in plot_appearance.lines:
        if line.vertical:
            ax.axvline(line.intercept, color=line.color,
                       linestyle=':', linewidth=2)
        else:
            ax.axhline(line.intercept, color=line.color,
                       linestyle=':', linewidth=2)

    # Activate the colorbar only if the Z axis contains data to plot
    if not isinstance(plot_data.z_data, EmptyAxisData):
        cbar = pyplot.colorbar(scatter, drawedges=False)
        cbar.ax.set_ylabel(plot_data.z_data.label, fontweight='bold')

    # Adjust layout
    pyplot.tight_layout()

    # Display NBD logo
    if not plot_appearance.hide_logo:
        import matplotlib.cbook as cbook
        import matplotlib.image as image
        from pele_platform.Utilities.Helpers import get_data_file_path

        with cbook.get_sample_data(get_data_file_path('NBD.png')) as file:
            im = image.imread(file)

        newax = fig.add_axes([0.79, 0.80, 0.2, 0.2], anchor='NE', zorder=1)
        newax.imshow(im)
        newax.axis('off')

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

    # Plot it
    pyplot.show()


def scatter_plot(pele_data, plot_data, plot_appearance, save_to):
    """
    It generates an scatter plot.

    Parameters
    ----------
    pele_data : a Pandas.DataFrame object
        The DataFrame that contains the PELE simulation data
    plot_data : a PlotData object
        The PlotData containing the information to plot in each
        axes
    plot_appearance : a PlotAppearance object
        The appearance settings for the plot
    save_to : str
        If it is not None, the resulting plot will be saved to this path
    """
    from matplotlib import pyplot

    # Set seaborn style
    pyplot.style.use("seaborn-ticks")

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

        cmap = get_colormap(plot_appearance.colormap_name)
        norm = pyplot.Normalize(z_min, z_max)
        color = None
    else:
        cmap = None
        norm = None
        color = COLORS[plot_appearance.plot_color][1]

    fig, ax = pyplot.subplots()

    scatter = ax.scatter(x_values, y_values, c=z_values, s=20,
                         cmap=cmap, norm=norm, color=color)

    if z_values is not None:
        cbar = pyplot.colorbar(scatter)
        cbar.ax.set_ylabel(plot_data.z_data.label, fontweight='bold')

    # Axes settings
    title = plot_appearance.title
    if title is None:
        title = "PELE Scatter Plot"
    ax.set_title(title, fontweight='bold', pad=20)
    ax.margins(0.05)
    ax.set_facecolor(plot_appearance.background_color)
    pyplot.ylabel(plot_data.y_data.label, fontweight='bold')
    pyplot.xlabel(plot_data.x_data.label, fontweight='bold')
    pyplot.xlim([plot_data.x_data.lowest, plot_data.x_data.highest])
    pyplot.ylim([plot_data.y_data.lowest, plot_data.y_data.highest])

    # Add lines
    for line in plot_appearance.lines:
        if line.vertical:
            ax.axvline(line.intercept, color=line.color,
                       linestyle=':', linewidth=2)
        else:
            ax.axhline(line.intercept, color=line.color,
                       linestyle=':', linewidth=2)

    # Adjust layout
    pyplot.tight_layout()

    # Display NBD logo
    if not plot_appearance.hide_logo:
        import matplotlib.cbook as cbook
        import matplotlib.image as image
        from pele_platform.Utilities.Helpers import get_data_file_path

        with cbook.get_sample_data(get_data_file_path('NBD.png')) as file:
            im = image.imread(file)

        newax = fig.add_axes([0.79, 0.80, 0.2, 0.2], anchor='NE', zorder=1)
        newax.imshow(im)
        newax.axis('off')

    # Plot it or save it
    if save_to is None:
        pyplot.show()
    else:
        pyplot.savefig(save_to, dpi=200)


def density_plot(pele_data, plot_data, plot_appearance, save_to):
    """
    It generates a density plot.

    Parameters
    ----------
    pele_data : a Pandas.DataFrame object
        The DataFrame that contains the PELE simulation data
    plot_data : a PlotData object
        The PlotData containing the information to plot in each
        axes
    plot_appearance : a PlotAppearance object
        The appearance settings for the plot
    save_to : str
        If it is not None, the resulting plot will be saved to this path
    """
    from matplotlib import pyplot
    import seaborn as sns

    sns.set_style("ticks")

    x_values = plot_data.get_xs_from_pele_data(pele_data).to_numpy()
    y_values = plot_data.get_ys_from_pele_data(pele_data).to_numpy()

    color1, color2 = COLORS[plot_appearance.plot_color]

    cmap = sns.dark_palette(color2, reverse=True, as_cmap=True)

    ax = sns.JointGrid(x=x_values, y=y_values)

    if plot_appearance.display_edges is False:
        markers_alpha = 0.7
    else:
        markers_alpha = 0.4
        ax.plot_joint(sns.kdeplot, cmap=cmap, shade=False,
                      n_levels=plot_appearance.n_levels)

    ax.plot_joint(sns.scatterplot, color=color1, edgecolor=color2,
                  marker='o', alpha=markers_alpha, s=20)

    sns.kdeplot(x=x_values, ax=ax.ax_marg_x, color=color1, shade=True,
                alpha=0.5, edgecolor=color2)

    sns.kdeplot(y=y_values, ax=ax.ax_marg_y, color=color1, shade=True,
                alpha=0.5, edgecolor=color2)

    # Axes settings
    title = plot_appearance.title
    if title is None:
        title = "PELE Density Plot"
    ax.fig.suptitle(title, fontweight='bold')
    ax.ax_joint.set_xlabel(plot_data.x_data.label, fontweight='bold')
    ax.ax_joint.set_ylabel(plot_data.y_data.label, fontweight='bold')
    ax.ax_marg_x.set_xlim(plot_data.x_data.lowest,
                          plot_data.x_data.highest)
    ax.ax_marg_y.set_ylim(plot_data.y_data.lowest,
                          plot_data.y_data.highest)

    # Add lines
    for line in plot_appearance.lines:
        if line.vertical:
            ax.ax_joint.axvline(line.intercept,
                                color=line.color,
                                linestyle=':', linewidth=2)
        else:
            ax.ax_joint.axhline(line.intercept,
                                color=line.color,
                                linestyle=':', linewidth=2)

    # Adjust layout
    pyplot.tight_layout()

    # Display NBD logo
    if not plot_appearance.hide_logo:
        import matplotlib.cbook as cbook
        import matplotlib.image as image
        from pele_platform.Utilities.Helpers import get_data_file_path

        with cbook.get_sample_data(get_data_file_path('NBD.png')) as file:
            im = image.imread(file)

        newax = ax.fig.add_axes([0.79, 0.80, 0.2, 0.2], anchor='NE', zorder=1)
        newax.imshow(im)
        newax.axis('off')

    # Plot it or save it
    if save_to is None:
        pyplot.show()
    else:
        pyplot.savefig(save_to, dpi=200)


def histogram_plot(pele_data, plot_data, plot_appearance, save_to):
    """
    It generates a histogram.

    Parameters
    ----------
    pele_data : a Pandas.DataFrame object
        The DataFrame that contains the PELE simulation data
    plot_data : a PlotData object
        The PlotData containing the information to plot in each
        axes
    plot_appearance : a PlotAppearance object
        The appearance settings for the plot
    save_to : str
        If it is not None, the resulting plot will be saved to this path
    """
    from matplotlib import pyplot
    import seaborn as sns

    sns.set_style("ticks")

    x_values = plot_data.get_xs_from_pele_data(pele_data).to_numpy()
    y_values = plot_data.get_ys_from_pele_data(pele_data).to_numpy()

    color1, color2 = COLORS[plot_appearance.plot_color]

    cmap = sns.dark_palette(color2, reverse=True, as_cmap=True)

    ax = sns.JointGrid(x=x_values, y=y_values)

    if plot_appearance.display_edges is False:
        markers_alpha = 0.7
    else:
        markers_alpha = 0.4
        ax.plot_joint(sns.kdeplot, cmap=cmap, shade=False,
                      n_levels=plot_appearance.n_levels)

    ax.plot_joint(sns.scatterplot, color=color1, edgecolor=color2,
                  marker='o', alpha=markers_alpha, s=20)

    sns.histplot(x=x_values, ax=ax.ax_marg_x, color=color1,
                 bins=plot_appearance.n_bins[0], stat='count',
                 alpha=0.5, edgecolor=color2)

    sns.histplot(y=y_values, ax=ax.ax_marg_y, color=color1,
                 bins=plot_appearance.n_bins[1], stat='count',
                 alpha=0.5, edgecolor=color2)

    # Axes settings
    title = plot_appearance.title
    if title is None:
        title = "PELE Histogram"
    ax.fig.suptitle(title, fontweight='bold')
    ax.ax_joint.set_xlabel(plot_data.x_data.label, fontweight='bold')
    ax.ax_joint.set_ylabel(plot_data.y_data.label, fontweight='bold')
    ax.ax_marg_x.set_xlim(plot_data.x_data.lowest,
                          plot_data.x_data.highest)
    ax.ax_marg_y.set_ylim(plot_data.y_data.lowest,
                          plot_data.y_data.highest)

    # Add lines
    for line in plot_appearance.lines:
        if line.vertical:
            ax.ax_joint.axvline(line.intercept,
                                color=line.color,
                                linestyle=':', linewidth=2)
        else:
            ax.ax_joint.axhline(line.intercept,
                                color=line.color,
                                linestyle=':', linewidth=2)

    # Adjust layout
    pyplot.tight_layout()

    # Display NBD logo
    if not plot_appearance.hide_logo:
        import matplotlib.cbook as cbook
        import matplotlib.image as image
        from pele_platform.Utilities.Helpers import get_data_file_path

        with cbook.get_sample_data(get_data_file_path('NBD.png')) as file:
            im = image.imread(file)

        newax = ax.fig.add_axes([0.79, 0.80, 0.2, 0.2], anchor='NE', zorder=1)
        newax.imshow(im)
        newax.axis('off')

    # Plot it or save it
    if save_to is None:
        pyplot.show()
    else:
        pyplot.savefig(save_to, dpi=200)


# Main workflow to be executed
if __name__ == "__main__":
    # Import external libraries
    import pandas as pd

    # Avoid Python2
    from pele_platform.Checker.python_version import check_python_version

    check_python_version()

    # Parse command-line arguments
    csv_file, results_folder, output_folder, report_name, \
        trajectory_name, plot_type, xdata, ydata, zdata, \
        xlowest, xhighest, ylowest, yhighest, zlowest, zhighest, \
        colormap, color,  with_edges, n_levels, n_bins, \
        background_color, vertical_lines, horizontal_lines, \
        filters, title, hide_logo, save_to = parse_args()

    # Print header
    from pele_platform.constants import constants
    print(constants.plotter_version_header.format(__version__))

    # Print parameters
    print_parameters(csv_file, results_folder, output_folder,
                     report_name, trajectory_name, plot_type,
                     xdata, ydata, zdata, xlowest, xhighest,
                     ylowest, yhighest, zlowest, zhighest,
                     colormap, color, with_edges, n_levels,
                     n_bins, background_color, vertical_lines,
                     horizontal_lines, filters, title, hide_logo,
                     save_to)

    # Get PELE data
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

    # Parse PELE data
    x_data = parse_axis_data(xdata, xlowest, xhighest, pele_data)
    y_data = parse_axis_data(ydata, ylowest, yhighest, pele_data)
    z_data = parse_axis_data(zdata, zlowest, zhighest, pele_data)

    # Ask axis data to the user if it is missing
    requested = False
    if isinstance(x_data, EmptyAxisData):
        x_data = request_axis_data('X', 'y', pele_data)
        requested = True
    if isinstance(y_data, EmptyAxisData):
        y_data = request_axis_data('Y', 'y', pele_data)
        requested = True
    if (isinstance(z_data, EmptyAxisData) and requested and
            plot_type.lower() != 'density'):
        z_data = request_axis_data('Z', 'n', pele_data, optional=True)

    for i, filter in enumerate(filters, start=1):
        print(f'-> Applying filter {i}')
        pele_data = filter.apply(pele_data)

    # Initialize plot data
    plot_data = PlotData(x_data, y_data, z_data)

    # Initialize plot appearance class
    plot_appearance = PlotAppearance(colormap_name=colormap,
                                     plot_color=color,
                                     background_color=background_color,
                                     display_edges=with_edges,
                                     n_levels=n_levels,
                                     n_bins=n_bins,
                                     lines=vertical_lines + horizontal_lines,
                                     title=title,
                                     hide_logo=hide_logo)

    # See if plot data can be plotted
    if not plot_data.is_plottable():
        raise ValueError('Aborted: not enough data to generate the plot.')

    # Generate the right plot type
    if save_to is None:
        print('-> Displaying plot. Close the window to exit.')
    else:
        print(f'-> Saving plot to {save_to}')

    if plot_type.lower() == 'interactive':
        interactive_plot(pele_data, plot_data, plot_appearance, save_to)
    elif plot_type.lower() == 'scatter':
        scatter_plot(pele_data, plot_data, plot_appearance, save_to)
    elif plot_type.lower() == 'density':
        density_plot(pele_data, plot_data, plot_appearance, save_to)
    elif plot_type.lower() == 'histogram':
        histogram_plot(pele_data, plot_data, plot_appearance, save_to)
