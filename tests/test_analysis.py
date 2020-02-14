import os
import glob
import shutil
import pele_platform.constants.constants as cs
import pele_platform.Analysis.plots as pt

simulation_path = os.path.join(os.path.dirname(cs.DIR), "tests/data/output")
REPORT_NAME = "report_"
TRAJ_NAME = "trajectory_"


def test_plot_two_metrics(simulation_path=simulation_path, report_name=REPORT_NAME, traj_name=TRAJ_NAME):
    output_folder="tmp/tests"
    if os.path.exists(output_folder):
        shutil.rmtree(output_folder)
    analysis = pt.PostProcessor(report_name, traj_name, simulation_path, 1)
    output = analysis.plot_two_metrics(5, 4, 3, output_folder=output_folder)
    assert os.path.exists(output)
    output = analysis.plot_two_metrics(5, 4, output_folder=output_folder)
    assert os.path.exists(output)
    shutil.rmtree("tmp")

    




def test_best_structs(simulation_path=simulation_path, report_name=REPORT_NAME, traj_name=TRAJ_NAME, n_structs=1):
    output_folder="tmp/tests/BestStructs"
    if os.path.exists(output_folder):
        shutil.rmtree(output_folder)
    analysis = pt.PostProcessor(report_name, traj_name, simulation_path, 1)
    analysis.top_poses(5, n_structs, output_folder)
    files = glob.glob(os.path.join(output_folder, "*"))
    assert len(files) == 1
    shutil.rmtree("tmp")



