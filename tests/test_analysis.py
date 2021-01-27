import os
import logging
import glob
import pytest
import shutil
import pele_platform.constants.constants as cs
import pele_platform.main as main
import pele_platform.Analysis.plots as pt

test_path = os.path.join(cs.DIR, "Examples")
simulation_path = "../pele_platform/Examples/analysis/data/output"
data = "data"
REPORT_NAME = "report_"
TRAJ_NAME = "trajectory_"
ANALYSIS_ARGS = os.path.join(test_path, "analysis/input.yaml")
ANALYSIS_FLAGS0 = os.path.join(test_path, "analysis/input_flags0.yaml")
ANALYSIS_FLAGS = os.path.join(test_path, "analysis/input_flags.yaml")
ANALYSIS_MAE_ARGS = os.path.join(test_path, "analysis/input_mae.yaml")


def test_plot_two_metrics(
    simulation_path=simulation_path, report_name=REPORT_NAME, traj_name=TRAJ_NAME
):
    output_folder = "tmp/plots"
    if os.path.exists(output_folder):
        shutil.rmtree(output_folder)
    analysis = pt.PostProcessor(report_name, traj_name, simulation_path, 1)
    output = analysis.plot_two_metrics(4, 5, 6, output_folder=output_folder)
    assert os.path.exists(output)
    output = analysis.plot_two_metrics(5, 6, output_folder=output_folder)
    assert os.path.exists(output)


def test_best_structs(
    simulation_path=simulation_path,
    report_name=REPORT_NAME,
    traj_name=TRAJ_NAME,
    n_structs=1,
):
    """
    Tests the retrieval of best_structures in PostProcessor.
    """
    output_folder = "tmp/tests/BestStructs"
    if os.path.exists(output_folder):
        shutil.rmtree(output_folder)
    analysis = pt.PostProcessor(report_name, traj_name, simulation_path, 1)
    analysis.logger = logging.getLogger("logger")
    analysis.top_poses(5, n_structs, output_folder)
    files = glob.glob(os.path.join(output_folder, "*"))
    assert len(files) == 1


def test_csv_move_folder(
    simulation_path=simulation_path,
    report_name=REPORT_NAME,
    traj_name=TRAJ_NAME,
    n_structs=1,
):
    output_folder = "copy_folder"
    shutil.copytree(simulation_path, output_folder)
    analysis = pt.PostProcessor(report_name, traj_name, output_folder, 1)
    analysis.logger = logging.getLogger("logger")
    analysis.top_poses(5, n_structs, output_folder)
    shutil.rmtree(output_folder)


def test_analysis_0flag(ext_args=ANALYSIS_FLAGS0):
    job = main.run_platform(ext_args)
    assert os.path.exists(
        "../pele_platform/Examples/analysis/data/results/Plots/numberOfAcceptedPeleSteps_currentEnergy_distance0_plot.png"
    )
    assert os.path.exists(
        "../pele_platform/Examples/analysis/data/results/Plots/numberOfAcceptedPeleSteps_currentEnergy_sasaLig_plot.png"
    )
    assert os.path.exists(
        "../pele_platform/Examples/analysis/data/results/Plots/distance0_currentEnergy_plot.png"
    )
    assert os.path.exists(
        "../pele_platform/Examples/analysis/data/results/Plots/sasaLig_currentEnergy_plot.png"
    )
    assert job.analysis_nclust == 1
    assert (
        len(glob.glob("../pele_platform/Examples/analysis/data/results/Plots/*.png"))
        == 4
    )


def test_analysis_flags(ext_args=ANALYSIS_FLAGS):
    main.run_platform(ext_args)
    assert os.path.exists(
        "../pele_platform/Examples/analysis/data/results/Plots/numberOfAcceptedPeleSteps_currentEnergy_distance0_plot.png"
    )
    assert os.path.exists(
        "../pele_platform/Examples/analysis/data/results/Plots/distance0_currentEnergy_plot.png"
    )
    assert (
        len(glob.glob("../pele_platform/Examples/analysis/data/results/Plots/*.png"))
        == 2
    )


@pytest.mark.parametrize("yaml", [ANALYSIS_ARGS, ANALYSIS_MAE_ARGS])
def test_analysis(yaml):
    """
    Full run of the platform with only_analysis flag to catch any errors. The second input extracts best structures
    and clusters as MAE files with the metrics as properties.
    """
    main.run_platform(yaml)


def test_cluster():
    output_folder = "clusters"
    n_clusts = 2
    if os.path.exists(output_folder):
        shutil.rmtree(output_folder)
    analysis = pt.PostProcessor(
        "report_",
        "trajectory_",
        "../pele_platform/Examples/analysis/data/xtc",
        2,
        topology="../pele_platform/Examples/analysis/data/xtc/topologies/topology_0.pdb",
        residue="L01",
    )
    analysis.logger = logging.getLogger("logger")
    analysis.retrive_data()
    analysis.cluster_poses(2, 5, output_folder, nclusts=n_clusts)
    assert os.path.exists(output_folder)
    assert len(glob.glob(os.path.join(output_folder, "clust*.pdb"))) == n_clusts - 1
