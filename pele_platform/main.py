# -*- coding: utf-8 -*-
"""
This is the main module designed to run the NBD platform from command line.
"""
from pele_platform.Errors import custom_errors
from pele_platform.Utilities.Helpers.launcher import Launcher
from pele_platform.context import context

__author__ = "Nostrum Biodiscovery"
__email__ = "pelesupport@nostrumbiodiscovery.com"
__license__ = "Apache-2.0"


def parse_args(args=None):
    """
    Command line parser.

    Parameters
    ----------
    args : List[str]
        The list of strings to parse. Default is an empty list. If empty,
        input arguments are retrieved from the command-line

    Returns
    -------
    input_yaml : str
        The path pointing to the input yaml file.
    """
    if args is None:  # to avoid mutable default argument
        args = []

    from argparse import ArgumentParser

    parser = ArgumentParser(description='Automatic platform to launch PELE simulations.')
    parser.add_argument('input_file', type=str, help='YAML input file.')
    parsed_args = parser.parse_args(args) if args else parser.parse_args()
    return parsed_args.input_file


def run_platform_from_yaml(input_yaml):
    """
    High level function to run the PELE platform from an input yaml file.

    It will:
    1) Parse the input.yaml
    2) Launch job
    3) Return job parameters

    Parameters
    ----------
    input_yaml : str
        The path pointing to the input yaml file

    Returns
    -------
    job_params : a Parameters object
        The corresponding Parameters object with the parameters of the
        simulation
    """
    from pele_platform.Utilities.Helpers import yaml_parser
    context.reset_parameters()
    context.yaml_parser = yaml_parser.YamlParser(input_yaml)

    try:
        context.yaml_parser.read()
    except AttributeError:
        raise custom_errors.WrongYamlFile(
            "Input file: {}".format(input_yaml)
            + " does not look like a correct yaml file")

    launcher = Launcher()
    job_parameters = launcher.launch_pele()

    return job_parameters


# Main workflow to be executed
if __name__ == "__main__":
    # Avoid backend issues in matplotlib
    import matplotlib
    matplotlib.use("Agg")

    # Avoid Python2
    from .Checker.python_version import check_python_version
    check_python_version()

    # Parse yaml file from command-line arguments
    yaml_file = parse_args()

    # Call platform runner
    job_params = run_platform_from_yaml(yaml_file)
