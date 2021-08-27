import glob
import os
import pytest
import shutil

from pele_platform.constants import constants as cs
from pele_platform.Utilities.Helpers import helpers
from pele_platform import main
from .utils import check_file

test_path = os.path.join(cs.DIR, "Examples")


LOCAL_ADAPTIVE = [
    '"type" : "inverselyProportional",',
    '"epsilon": 0.25,',
    '"iterations" : 10,',
    '"peleSteps" : 50,',
    '"values" : [2.0, 4, 6],',
    '"conditions": [1, 0.5, 0.0]',
]

GLOBAL_ADAPTIVE = [
    '"epsilon": 0.25,',
    '"type" : "inverselyProportional",',
    '"iterations" : 50,',
    '"peleSteps" : 12,',
    '"values" : [2.5, 4, 6],',
    '"conditions": [1, 0.5, 0.0]',
]

PELE = [
    '"sideChainPredictionFrequency" : 2,',
    '"temperature": 1500,',
    '"numberOfStericTrials": 200,',
    '"overlapFactor": 0.65',
]


def test_site_finder_skipref():
    """
    Tests skip_refinement flag, so that site_finder runs global exploration only.
    """
    yaml_file = os.path.join(test_path, "site_finder/input_skipref.yaml")
    job, = main.run_platform_from_yaml(yaml_file)
    refinement_simulation = os.path.join(
        os.path.dirname(job.pele_dir), "Clusters"
    )
    assert not os.path.exists(refinement_simulation)


@pytest.mark.parametrize("yaml", ["site_finder/input_global.yaml", "site_finder/input_global_xtc.yaml"])
def test_site_finder_production(yaml):
    """
    Test site finder production on both XTC and PDB trajectories.
    """
    yaml_file = os.path.join(test_path, yaml)
    job = main.run_platform_from_yaml(yaml_file)

    randomization_output = glob.glob(os.path.join(job[0].inputs_dir, "input*pdb"))

    refinement_input = glob.glob(
        os.path.join(os.path.dirname(job[0].pele_dir), "2_Clusters/cluster_A.pdb")
    )

    files_refinement_output = glob.glob(os.path.join(job[2].pele_dir, "results/top_poses/*.pdb"))

    best_energy_input = os.path.join(
        os.path.dirname(job[0].pele_dir),
        "2_Clusters",
        "cluster_A.pdb"
    )

    assert len(randomization_output) == job[0].poses
    assert best_energy_input in refinement_input
    assert files_refinement_output


def test_working_folder(output="site_finder"):
    """
    Tests custom working folder.
    """
    yaml_file = os.path.join(test_path, "site_finder/input_folder.yaml")
    helpers.check_remove_folder(output)
    job, = main.run_platform_from_yaml(yaml_file)
    assert os.path.exists(job.folder)
    helpers.check_remove_folder(output)


def test_site_finder_restart():
    """
    Tests if site finder can pick up adaptive_restart flag and sort out directories correctly when it's not a real
    adaptive restart (output folder does not exist yet).
    """
    expected_path = "site_finder_adaptive_restart"
    helpers.check_remove_folder(expected_path)
    shutil.copytree(os.path.join(test_path, "site_finder", "site_finder"), expected_path)

    yaml_file = os.path.join(test_path, "site_finder/input_restart.yaml")
    job = main.run_platform_from_yaml(yaml_file)

    files_refinement_input = glob.glob(os.path.join(os.path.dirname(job[0].pele_dir), "2_Clusters/*.pdb"))
    files_refinement_output = glob.glob(os.path.join(job[2].pele_dir, "results/top_poses/*.pdb"))

    assert files_refinement_input
    assert files_refinement_output

    helpers.check_remove_folder(expected_path)


@pytest.mark.parametrize(("yaml_file", "adaptive_lines"),
                         [
                             ("local.yaml", LOCAL_ADAPTIVE),
                             ("global.yaml", GLOBAL_ADAPTIVE)
                         ])
def test_defaults(yaml_file, adaptive_lines):
    """
    Tests default adaptive and pele parameters for both global and local site finder exploration.
    """
    yaml_file = os.path.join(test_path, "site_finder", yaml_file)
    params, = main.run_platform_from_yaml(yaml_file)

    errors = check_file(params.pele_dir, "pele.conf", PELE, errors=[])
    errors = check_file(params.pele_dir, "adaptive.conf", adaptive_lines, errors)
    assert not errors
