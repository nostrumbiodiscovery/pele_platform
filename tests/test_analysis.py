import os
import glob
import pytest
import shutil
import pele_platform.constants.constants as cs
import pele_platform.main as main
from pele_platform.analysis import Analysis, DataHandler, Plotter

test_path = os.path.join(cs.DIR, "Examples")
simulation_path = "../pele_platform/Examples/analysis/data/output"
data = "data"
REPORT_NAME = "report"
TRAJ_NAME = "trajectory.pdb"
ANALYSIS_ARGS = os.path.join(test_path, "analysis/input.yaml")
ANALYSIS_FLAGS0 = os.path.join(test_path, "analysis/input_flags0.yaml")
ANALYSIS_FLAGS = os.path.join(test_path, "analysis/input_flags.yaml")
ANALYSIS_XTC_ARGS = os.path.join(test_path, "analysis/input_xtc.yaml")


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
    if os.path.exists(output_folder):
        shutil.rmtree(output_folder)

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
def test_top_poses(n_poses, expected_energies):
    """
    Checks if data_handler extracts the correct number of top poses and associated metrics.
    Returns
    -------
        Folder with top poses.
    """
    output_folder = "tmp/top_poses"
    if os.path.exists(output_folder):
        shutil.rmtree(output_folder)

    analysis = Analysis(
        resname="LIG",
        chain="Z",
        simulation_output="../pele_platform/Examples/clustering",
        skip_initial_structures=False,
    )
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
    Runs full simulation with input.yaml with some unusual flags, check the number of created plots and their names to
    ensure correct metrics were take into account.
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
        Folder with plots
    """
    output_folder = "../pele_platform/Examples/analysis/data/results/plots/"
    if os.path.exists(output_folder):
        shutil.rmtree(output_folder)

    main.run_platform(yaml_file)

    for file in expected_files:
        file_path = os.path.join(output_folder, file)
        assert os.path.exists(file_path)

    all_files = glob.glob(os.path.join(output_folder, "*png"))
    assert len(all_files) == n_expected_outputs


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
    shutil.rmtree(results_folder)


@pytest.mark.parametrize(
    ("method", "bandwidth", "n_clusters"),
    [
        ("hdbscan", 5, 0),  # only gets orphan clusters [-1]
        ("meanshift", 100, 1),
        ("meanshift", 30, 3),
        ("gaussianmixture", 1, 2),
    ],
)
def test_clustering_methods(method, bandwidth, n_clusters):
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

    if os.path.exists(working_folder):
        shutil.rmtree(working_folder)

    analysis = Analysis(
        resname="LIG",
        chain="Z",
        simulation_output="../pele_platform/Examples/clustering",
        skip_initial_structures=False,
        bandwidth=bandwidth,
        analysis_nclust=n_clusters,
        clustering_method=method,
    )
    analysis.generate_clusters(working_folder, method)
    assert len(glob.glob(results)) == n_clusters


def test_analysis_api():
    """
    Runs full analysis workflow (with GMM clustering).
    Returns
    -------
        Returns a directory with top_poses, clusters and plots.
    """
    working_folder = "full_analysis"
    output = "../pele_platform/Examples/clustering"
    n_clusts = 3

    analysis = Analysis(
        resname="LIG",
        chain="Z",
        simulation_output=output,
        skip_initial_structures=False,
        analysis_nclust=n_clusts,
    )
    analysis.generate(working_folder, "gaussianmixture")

    # Check if reports exist
    assert os.path.exists(os.path.join(working_folder, "data.csv"))
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

    # Check if data.csv exists and is not empty
    data_csv = os.path.join(working_folder, "data.csv")
    assert os.path.exists(data_csv)

    with open(data_csv, "r") as file:
        lines = file.readlines()
        assert len(lines) == 8
        assert lines[0] == "Step,numberOfAcceptedPeleSteps,currentEnergy,Binding Energy,sasaLig,epoch,trajectory,Cluster\n"
