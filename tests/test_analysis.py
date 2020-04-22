import os
import glob
import shutil
import pele_platform.constants.constants as cs
import pele_platform.main as main
import pele_platform.Analysis.plots as pt

test_path = os.path.join(cs.DIR, "Examples")
simulation_path = "data/output"
REPORT_NAME = "report_"
TRAJ_NAME = "trajectory_"
ANALYSIS_ARGS = os.path.join(test_path, "analysis/input.yaml")
ANALYSIS_FLAGS0 = os.path.join(test_path, "analysis/input_flags0.yaml")
ANALYSIS_FLAGS = os.path.join(test_path, "analysis/input_flags.yaml")
ANALYSIS_MAE_ARGS = os.path.join(test_path, "analysis/input_mae.yaml")


def test_plot_two_metrics(simulation_path=simulation_path, report_name=REPORT_NAME, traj_name=TRAJ_NAME):
    output_folder="tmp/plots"
    if os.path.exists(output_folder):
        shutil.rmtree(output_folder)
    analysis = pt.PostProcessor(report_name, traj_name, simulation_path, 1)
    output = analysis.plot_two_metrics(4, 5, 6, output_folder=output_folder)
    assert os.path.exists(output)
    output = analysis.plot_two_metrics(5, 6, output_folder=output_folder)
    assert os.path.exists(output)

def test_best_structs(simulation_path=simulation_path, report_name=REPORT_NAME, traj_name=TRAJ_NAME, n_structs=1):
    output_folder="tmp/tests/BestStructs"
    if os.path.exists(output_folder):
        shutil.rmtree(output_folder)
    analysis = pt.PostProcessor(report_name, traj_name, simulation_path, 1)
    analysis.top_poses(5, n_structs, output_folder)
    files = glob.glob(os.path.join(output_folder, "*"))
    assert len(files) == 1

def test_analysis_0flag(ext_args=ANALYSIS_FLAGS0):
    main.run_platform(ext_args)
    assert os.path.exists("data/results/Plots/numberOfAcceptedPeleSteps_currentEnergy_distance0_plot.png")
    assert os.path.exists("data/results/Plots/numberOfAcceptedPeleSteps_currentEnergy_sasaLig_plot.png")
    assert os.path.exists("data/results/Plots/distance0_currentEnergy_plot.png")
    assert os.path.exists("data/results/Plots/sasaLig_currentEnergy_plot.png")
    assert len(glob.glob("data/results/Plots/*.png")) == 4

def test_analysis_flag(ext_args=ANALYSIS_FLAGS):
    main.run_platform(ext_args)
    assert os.path.exists("data/results/Plots/numberOfAcceptedPeleSteps_currentEnergy_distance0_plot.png")
    assert os.path.exists("data/results/Plots/distance0_currentEnergy_plot.png")
    assert len(glob.glob("data/results/Plots/*.png")) == 2

def test_analysis(ext_args=ANALYSIS_ARGS):
    main.run_platform(ext_args)

def test_analysis_mae(ext_args=ANALYSIS_MAE_ARGS):
    os.system("rm data/*/*summary*")
    main.run_platform(ext_args)
