import os
import glob
import pytest

import pele_platform.constants.constants as cs
import pele_platform.main as main
from pele_platform.analysis import Analysis, DataHandler, Plotter
from pele_platform.Utilities.Helpers.helpers import check_remove_folder
from . import test_adaptive as ta

test_path = os.path.join(cs.DIR, "Examples")
simulation_path = "../pele_platform/Examples/analysis/data/output"
REPORT_NAME = "report"
TRAJ_NAME = "trajectory.pdb"

ANALYSIS_ARGS = os.path.join(test_path, "analysis/input.yaml")
ANALYSIS_FLAGS0 = os.path.join(test_path, "analysis/input_flags0.yaml")
ANALYSIS_FLAGS = os.path.join(test_path, "analysis/input_flags.yaml")
ANALYSIS_XTC_ARGS = os.path.join(test_path, "analysis/input_xtc.yaml")

expected_energies = [-8699.17, -8690.68, -8711.55, -8709.73, -8687.3, -8706.14, -8712.07]


@pytest.mark.parametrize(("x", "y", "z"), [(4, 5, 6), (5, 6, None)])
def test_plotter(x, y, z):
    """
    Checks if the scatter and KDE plots are created correctly.
    Parameters
    ----------
    x : int
        Metric to x
    y :
        Metric to y
    z :
        Metric to z
    Returns
    -------
        Folder with plots.
    """
    output_folder = "tmp/plots"
    check_remove_folder(output_folder)

    data_handler = DataHandler(
        sim_path=simulation_path,
        report_name=REPORT_NAME,
        trajectory_name=TRAJ_NAME,
        be_column=5,
    )
    dataframe = data_handler.get_reports_dataframe()
    plotter = Plotter(dataframe)
    output_scatter = plotter.plot_two_metrics(x, y, z, output_folder=output_folder)
    output_kde = plotter.plot_kde(x, y, output_folder=output_folder, kde_structs=10)

    assert os.path.exists(output_scatter)
    assert os.path.exists(output_kde)


@pytest.mark.parametrize(
    ("n_poses", "expected_energies"),
    [
        (0, []),
        (1, [0.879]),
        (4, [0.879, 2.203, 3.563, 6.624]),
    ],
)
def test_generate_top_poses(analysis, n_poses, expected_energies):
    """
    Checks if data_handler extracts the correct number of top poses and associated metrics.
    Returns
    -------
        Folder with top poses.
    """
    output_folder = "tmp/top_poses"
    check_remove_folder(output_folder)

    top_poses = analysis.generate_top_poses(output_folder, n_poses)
    top_poses_rounded = [round(pose, 3) for pose in top_poses]

    # Check if correct energy values were extracted
    assert len(top_poses) == n_poses
    for energy in expected_energies:
        assert energy in top_poses_rounded

    # Check if correct number of files was saved
    results = [
        os.path.basename(file)
        for file in glob.glob(os.path.join(output_folder, "*pdb"))
    ]
    assert len(results) == n_poses


@pytest.mark.parametrize(
    ("yaml_file", "n_expected_outputs", "expected_files"),
    [
        (
                ANALYSIS_FLAGS0,
                4,
                [
                    "distance0_Binding_Energy_plot.png",
                    "currentEnergy_Binding_Energy_distance0_plot.png",
                    "sasaLig_Binding_Energy_plot.png",
                    "currentEnergy_Binding_Energy_sasaLig_plot.png",
                ],
        ),
        (
                ANALYSIS_FLAGS,
                2,
                [
                    "currentEnergy_Binding_Energy_distance0_plot.png",
                    "distance0_Binding_Energy_plot.png",
                ],
        ),
    ],
)
def test_analysis_flags(yaml_file, n_expected_outputs, expected_files):
    """
    Runs full simulation with input.yaml with some unusual flags, check the number of top poses, created plots and their
    names to ensure correct metrics were take into account.
    Parameters
    ----------
    yaml_file : str
        Path to input.yaml
    n_expected_outputs : int
        Number of expected plots.
    expected_files : List[str]
        List of expected plot names.
    Returns
    -------
        Folder with results.
    """
    output_folder = "../pele_platform/Examples/analysis/data/results"
    plots_folder = os.path.join(output_folder, "plots")
    top_poses_folder = os.path.join(output_folder, "top_poses", "*pdb")

    check_remove_folder(output_folder)

    main.run_platform(yaml_file)

    # Check if all expected file names are present
    for file in expected_files:
        file_path = os.path.join(plots_folder, file)
        assert os.path.exists(file_path)

    # Check number of created plots and top poses
    all_plots = glob.glob(os.path.join(plots_folder, "*png"))
    assert len(all_plots) == n_expected_outputs

    all_top_poses = glob.glob(top_poses_folder)
    assert len(all_top_poses) == 0


@pytest.mark.parametrize(("yaml_file", "expected_poses", "expected_clusters"),
                         [(ANALYSIS_ARGS, 1, 0), (ANALYSIS_XTC_ARGS, 3, 2)])
def test_analysis_production(yaml_file, expected_poses, expected_clusters):
    """
    Runs production analysis from input.yaml, both for PDB and XTC trajectories.
    Parameters
    ----------
    yaml_file : str
        Path to input.yaml
    Returns
    -------
        Parameters object with simulation parameters.
    """
    job_params = main.run_platform(yaml_file)

    results_folder = os.path.join(job_params.pele_dir, "results")
    top_poses = glob.glob(os.path.join(results_folder, "top_poses/*pdb"))
    clusters = glob.glob(os.path.join(results_folder, "clusters/*pdb"))

    assert len(top_poses) == expected_poses
    assert len(clusters) == expected_clusters

    # Clean up
    check_remove_folder(results_folder)


@pytest.mark.parametrize(
    ("method", "bandwidth", "n_clusters"),
    [
        ("hdbscan", 5, 0),  # only gets orphan clusters [-1]
        ("meanshift", 100, 1),
        ("meanshift", 30, 3),
        ("gaussianmixture", 1, 2),
    ],
)
def test_generate_clusters(analysis, method, bandwidth, n_clusters):
    """
    Checks if built-in clustering methods are producing expected number of clusters.

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
        Folder with clusters, plots and report.
    """
    working_folder = "clustering_method"
    results = os.path.join(working_folder, "*pdb")
    check_remove_folder(working_folder)

    analysis.generate_clusters(working_folder, method, bandwidth=bandwidth, analysis_nclust=n_clusters)
    assert len(glob.glob(results)) == n_clusters


def test_api_analysis_generation(analysis):
    """
    Runs full analysis workflow (with GMM clustering).
    Returns
    -------
        Returns a directory with top_poses, clusters and plots.
    """
    working_folder = "full_analysis"
    check_remove_folder(working_folder)
    n_clusts = 3
    analysis.generate(working_folder, "gaussianmixture", analysis_nclust=n_clusts)

    # Check if reports exist
    assert os.path.exists(os.path.join(working_folder, "summary.pdf"))

    # Check plots
    plots = glob.glob(os.path.join(working_folder, "plots", "*png"))
    assert len(plots) == 2

    # Check top poses
    top_poses = glob.glob(os.path.join(working_folder, "top_poses", "*pdb"))
    assert len(top_poses) == 7

    # Check clusters
    clusters = glob.glob(os.path.join(working_folder, "clusters", "*pdb"))
    assert len(clusters) == n_clusts

    # Check cluster representatives CSV by testing for the presence of a representative line
    errors = ta.check_file(os.path.join(working_folder, "clusters"), "top_selections.csv",
                           "1,B,../pele_platform/Examples/clustering/0/trajectory_3.pdb,0", [])
    assert not errors

    # Check if data.csv exists and is not empty
    data_csv = os.path.join(working_folder, "data.csv")
    assert os.path.exists(data_csv)

    with open(data_csv, "r") as file:
        lines = file.readlines()
        assert len(lines) == 8
        assert lines[0] == "Step,numberOfAcceptedPeleSteps,currentEnergy,Binding Energy,sasaLig,epoch,trajectory," \
                           "Cluster\n"


def test_check_existing_directory(generate_folders):
    """
    Checks if tester of existing dir
    Parameters
    ----------
    generate_folders : pytest.fixture
        Pytest fixture that generates "results" folders for testing.
    Returns
    -------
        A new "results" folder.
    """
    new_path = Analysis._check_existing_directory("results")
    assert new_path == "results_3"

    folders = glob.glob("results*")
    check_remove_folder(*folders)


@pytest.mark.parametrize("max_coordinates", [1, 3])
def test_extract_and_filter_coordinates(analysis, max_coordinates):
    """
    Checks coordinates and dataframe extraction and its subsequent filtering.
    Parameters
    ----------
    analysis : Analysis
        Analysis object created in a fixture.
    max_coordinates : int
        Number of coordinates to extract per ligand.
    Returns
    -------

    """
    coordinates, dataframe = analysis._extract_coordinates(max_coordinates)
    coordinates_filtered, _, _ = analysis._filter_coordinates(coordinates, dataframe, 0.5)

    assert len(coordinates) == 7
    assert coordinates.shape == (7, max_coordinates, 3)  # 7 poses, n atoms, 3 coordinates

    assert len(dataframe) == 7
    assert dataframe["currentEnergy"].tolist().sort() == expected_energies.sort()
    assert len(coordinates_filtered) == 4


def test_extract_poses(analysis):
    """
    Tests poses extraction from dataframe.
    Parameters
    ----------
    analysis : Analysis object
        Created in analysis fixture.
    Returns
    -------
        Energy values for extracted poses and a folder with PDB files.
    """

    output = "extracted_poses"
    check_remove_folder(output)

    values = analysis._extract_poses(analysis._dataframe, "currentEnergy", output)
    poses = glob.glob(os.path.join(output, "*pdb"))

    assert values.sort() == expected_energies.sort()
    assert len(poses) == 7


@pytest.mark.parametrize(("filter", "threshold", "expected_length"),
                         [(True, 0.5, 1),
                          (False, None, 7)])
def test_get_dataframe(analysis, filter, threshold, expected_length):
    """
    Tests dataframe filtering based on highest energy. Filtering happens twice, for binding energy and total energy,
    therefore a threshold of 50% will only return one pose out of 7.
    Parameters
    ----------
    analysis : Analysis
        Analysis object created in a fixture.
    Returns
    -------
        Filtered dataframe.
    """

    df = analysis.get_dataframe(filter=filter, threshold=threshold)
    assert len(df) == expected_length


@pytest.fixture
def generate_folders():
    """
    Generates a few random folders to test _check_existing_directory.
    Returns
    -------
        Three "results" folders.
    """

    folder_template = "results_{}"
    folders = [folder_template.format(i) for i in range(3)] + ["results"]

    for folder in folders:
        os.mkdir(folder)


@pytest.fixture
def analysis():
    """
    Creates Analysis object to use in other tests.
    Returns
    -------
        Analysis object.
    """
    output = "../pele_platform/Examples/clustering"

    analysis = Analysis(resname="LIG", chain="Z",
                        simulation_output=output,
                        skip_initial_structures=False)

    return analysis
