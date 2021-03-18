import shutil
import glob
import os
import pytest

from pele_platform.constants import constants as cs
from pele_platform import main

test_path = os.path.join(cs.DIR, "Examples")


def test_site_finder_skipref():
    yaml_file = os.path.join(test_path, "site_finder/input_skipref.yaml")
    job, job2 = main.run_platform(yaml_file)
    refinement_simulation = os.path.join(
        os.path.dirname(job.pele_dir), "2_refinement_simulation"
    )

    assert not os.path.exists(refinement_simulation)
    assert not job2


def test_site_finder_pdb():
    yaml_file = os.path.join(test_path, "site_finder/input_global.yaml")
    job, job2 = main.run_platform(yaml_file)

    # checkpoints
    refinement_input = glob.glob(
        os.path.join(os.path.dirname(job.pele_dir), "refinement_input/*.pdb")
    )
    nfiles_refinement = len(
        glob.glob(os.path.join(job2.pele_dir, "results/top_poses/*.pdb"))
    )
    best_energy_input = os.path.join(
        os.path.dirname(job.pele_dir),
        "refinement_input",
        "epoch0_trajectory_3.1_BindingEnergy-1.24648.pdb",
    )

    # test
    assert best_energy_input in refinement_input
    assert nfiles_refinement > 0


def test_site_finder_xtc():
    yaml_file = os.path.join(test_path, "site_finder/input_global_xtc.yaml")
    job, job2 = main.run_platform(yaml_file)

    # checkpoints
    best_energy_input = os.path.join(
        os.path.dirname(job.pele_dir),
        "refinement_input",
        "epoch0_trajectory_3.1_BindingEnergy-1.24648.pdb",
    )
    refinement_input = glob.glob(
        os.path.join(os.path.dirname(job.pele_dir), "refinement_input/*.pdb")
    )
    nfiles_refinement = len(
        glob.glob(os.path.join(job2.pele_dir, "results/top_poses/*.pdb"))
    )

    # checking if all temporary top_poses (needed to select refinement input) were removed
    temp_bs_dir = os.path.join(job.pele_dir, job.output, "top_poses")

    # test
    assert best_energy_input in refinement_input
    assert nfiles_refinement > 0
    assert not os.path.exists(temp_bs_dir)


def test_working_folder(output="site_finder"):
    yaml_file = os.path.join(test_path, "site_finder/input_folder.yaml")
    if os.path.exists(output):
        shutil.rmtree(output)
    job, _ = main.run_platform(yaml_file)
    assert os.path.exists(job.folder)


@pytest.mark.skip(reason="Implemented in pele_platform 2.0.")
def test_site_finder_restart():
    yaml_file = os.path.join(test_path, "site_finder/input_restart.yaml")
    job, job2 = main.run_platform(yaml_file)

    nfiles = len(
        glob.glob(os.path.join(os.path.dirname(job.pele_dir), "refinement_input/*.pdb"))
    )
    nfiles_refinement = len(
        glob.glob(os.path.join(job2.pele_dir, "results/top_poses/*.pdb"))
    )
    # test
    assert nfiles == job.n_components
    assert nfiles_refinement
