import os
import logging
import glob
import pytest
import shutil
import pele_platform.constants.constants as cs
import pele_platform.main as main
import pele_platform.analysis.plot as pt
from pele_platform.Utilities.Parameters.parameters import ParametersBuilder
from pele_platform.analysis import DataHandler, Plotter
from pele_platform.analysis import (GaussianMixtureClustering,
                                    HDBSCANClustering,
                                    MeanShiftClustering)

test_path = os.path.join(cs.DIR, "Examples")
simulation_path = "../pele_platform/Examples/analysis/data/output"
data = "data"
REPORT_NAME = "report"
TRAJ_NAME = "trajectory.pdb"
ANALYSIS_ARGS = os.path.join(test_path, "analysis/input.yaml")
ANALYSIS_FLAGS0 = os.path.join(test_path, "analysis/input_flags0.yaml")
ANALYSIS_FLAGS = os.path.join(test_path, "analysis/input_flags.yaml")
ANALYSIS_XTC_ARGS = os.path.join(test_path, "analysis/input_xtc.yaml")


# @pytest.fixture
# def parameters():
#
#     builder = ParametersBuilder()
#     parameters = builder.build_adaptive_variables()


@pytest.mark.parametrize(("x", "y", "z"), [(4, 5, 6), (5, 6, None)])
def test_plot_two_metrics(x, y, z):

    output_folder = "tmp/plots"
    if os.path.exists(output_folder):
        shutil.rmtree(output_folder)

    data_handler = DataHandler(sim_path=simulation_path, report_name=REPORT_NAME, trajectory_name=TRAJ_NAME, be_column=5)
    dataframe = data_handler.get_reports_dataframe()
    plotter = Plotter(dataframe)
    output = plotter.plot_two_metrics(x, y, z, output_folder=output_folder)
    assert os.path.exists(output)


def test_best_structs():
    """
    Checks if data_handler correctly extracts metrics for best structures.
    Returns
    -------

    """
    metric = "Binding Energy"
    n_poses = 1
    output_folder = "tmp/tests/BestStructs"
    if os.path.exists(output_folder):
        shutil.rmtree(output_folder)

    data_handler = DataHandler(sim_path=simulation_path, report_name=REPORT_NAME, trajectory_name=TRAJ_NAME,
                               be_column=5)
    top_poses = data_handler.get_top_entries(metric, n_poses)
    assert len(top_poses[metric]) == n_poses


def test_csv_move_folder():
    # TODO: Is it still necessary?
    metric = "Binding Energy"
    n_poses = 1

    output_folder = "copy_folder"
    if os.path.exists(output_folder):
        shutil.rmtree(output_folder)

    shutil.copytree(simulation_path, output_folder)
    data_handler = DataHandler(sim_path=simulation_path, report_name=REPORT_NAME, trajectory_name=TRAJ_NAME,
                               be_column=5)
    top_poses = data_handler.get_top_entries(metric, n_poses)


def test_analysis_0flag(ext_args=ANALYSIS_FLAGS0):
    if os.path.exists("../pele_platform/Examples/analysis/data/results/plots/"):
        shutil.rmtree("../pele_platform/Examples/analysis/data/results/plots/")

    job = main.run_platform(ext_args)
    assert os.path.exists(
        "../pele_platform/Examples/analysis/data/results/plots/distance0_Binding_Energy_plot.png"
    )
    assert os.path.exists(
        "../pele_platform/Examples/analysis/data/results/plots/currentEnergy_Binding_Energy_distance0_plot.png"
    )
    assert os.path.exists(
        "../pele_platform/Examples/analysis/data/results/plots/sasaLig_Binding_Energy_plot.png"
    )
    assert os.path.exists(
        "../pele_platform/Examples/analysis/data/results/plots/currentEnergy_Binding_Energy_sasaLig_plot.png"
    )
    assert job.analysis_nclust == 4
    assert (
        len(glob.glob("../pele_platform/Examples/analysis/data/results/plots/*.png"))
        == 4
    )


def test_analysis_flag(ext_args=ANALYSIS_FLAGS):
    if os.path.exists("../pele_platform/Examples/analysis/data/results/plots/"):
        shutil.rmtree("../pele_platform/Examples/analysis/data/results/plots/")

    main.run_platform(ext_args)
    assert os.path.exists(
        "../pele_platform/Examples/analysis/data/results/plots/currentEnergy_Binding_Energy_distance0_plot.png"
    )
    assert os.path.exists(
        "../pele_platform/Examples/analysis/data/results/plots/distance0_Binding_Energy_plot.png"
    )
    assert (
        len(glob.glob("../pele_platform/Examples/analysis/data/results/plots/*.png"))
        == 2
    )


def test_analysis_production(ext_args=ANALYSIS_ARGS):
    main.run_platform(ext_args)


def test_analysis_xtc(ext_args=ANALYSIS_XTC_ARGS):
    from pele_platform.Utilities.Helpers import yaml_parser
    from pele_platform.Checker import valid_flags
    from pele_platform.Utilities.Helpers.launcher import Launcher

    yaml_obj = yaml_parser.YamlParser(ext_args,
                                      valid_flags.VALID_FLAGS_PLATFORM)

    yaml_obj.read()

    # Modify this value for this specific case
    yaml_obj.analysis_nclust = 3

    # Initialize job launcher
    launcher = Launcher(yaml_obj)

    # Run launcher
    job_params = launcher.launch()

    # Return job parameters
    return job_params

def test_cluster_default():
    n_clusts = 2

    data_handler = DataHandler(sim_path='../pele_platform/Examples/clustering',
                               report_name=REPORT_NAME,
                               trajectory_name=TRAJ_NAME,
                               be_column=5)

    coordinates, dataframe = data_handler.extract_raw_coords('LIG')

    clustering = GaussianMixtureClustering(2)
    labels = clustering.get_clusters(coordinates)

    assert len(set(labels)) == n_clusts


@pytest.mark.parametrize(("method", "bandwidth", "n_clusters"), [
    ("hdbscan", 20, 1),
    ("meanshift", 100, 1),
    ("meanshift", 10, 3)])
def test_clustering_methods(method, bandwidth, n_clusters):
    """
    Checks built-in clustering methods and report generation.

    Parameters
    ----------
    method : str
        Built-in clustering method, e.g. "dbscan".
    bandwidth : float
        Bandwidth for meanshift (or epsilon for DBSCAN).
    n_clusters : int
        Number of clusters for the Gaussian mixture model.

    Returns
    -------
        Folder with clusters and report.
    """

    data_handler = DataHandler(sim_path='../pele_platform/Examples/clustering',
                               report_name=REPORT_NAME,
                               trajectory_name=TRAJ_NAME,
                               be_column=5)

    coordinates, dataframe = data_handler.extract_raw_coords('LIG')

    if method == 'hdbscan':
        clustering = HDBSCANClustering(bandwidth)
    elif method == 'meanshift':
        clustering = MeanShiftClustering(bandwidth)
    else:
        assert False

    labels = clustering.get_clusters(coordinates)

    assert len(set(labels)) == n_clusters
