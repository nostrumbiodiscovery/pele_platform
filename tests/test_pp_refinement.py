import os
import pytest
import yaml

from pele_platform.constants import constants
from pele_platform import main
from .test_adaptive import check_file

test_path = os.path.join(constants.DIR, "Examples")

PELE = [
    '"anmFrequency" : 0,',
    '"sideChainPredictionFrequency" : 1,',
    '"minimizationFrequency" : 1,',
    '"sideChainPredictionRegionRadius" : 6,',
    '"perturbationCOMConstraintConstant" : 0,',
    '"activateProximityDetection": false,',
    '"perturbationType":"naive",',
    '"rotationAngles": "nonCoupled",',
    '"algorithm": "CARTESIANS", "nodes": { "atoms": { "names": [ "_CA_" ]} },',
    '"algorithm" : "zhexin",',
    '"algorithm" : "TruncatedNewton",',
]

ADAPTIVE = [
    '"type" : "epsilon",',
    '"epsilon": 0.2,',
    '"iterations" : 1,',
    '"peleSteps" : 250,',
    '"runEquilibration" : false,',
    '"equilibrationMode": "equilibrationSelect",',
    '"equilibrationLength" : 1,',
    '"values" : [3.0, 5, 6],',
    '"conditions": [0.2, 0.1, 0]',
    '"type" : "heaviside",',
    '"ligandChain": "A",',
    '"alternativeStructure" : true,',
    '"contactThresholdDistance" : 8',
]


def test_defaults():
    """
    Runs in debug mode to check all the package defaults in pele.conf and adaptive.conf.
    """
    yaml_file = os.path.join(test_path, "pp_refinement", "input_defaults.yaml")
    parameters = main.run_platform_from_yaml(yaml_file)

    errors = check_file(parameters.pele_dir, "pele.conf", PELE, [])
    errors = check_file(parameters.pele_dir, "adaptive.conf", ADAPTIVE, errors)
    assert not errors


@pytest.mark.parametrize(
    "pdb_file",
    [
        "../pele_platform/Examples/pp_refinement/4cc*.pdb",
        "../pele_platform/Examples/pp_refinement/6UD7_C_0.pdb",
        "../pele_platform/Examples/pp_refinement/DCAF1*.pdb",
        "../pele_platform/Examples/pp_refinement/light*.pdb",
        "../pele_platform/Examples/pp_refinement/model*.pdb",
    ],
)
def test_production(pdb_file):
    """
    Runs the whole workflow using several structures originating from different software by editing the YAML file
    between runs.

    Parameters
    ----------
    pdb_file : str
        Path to the PDB file with system.
    """
    yaml_file = os.path.join(test_path, "pp_refinement", "input.yaml")

    with open(yaml_file, "r") as file:
        user_dict = yaml.load(file)

    user_dict["system"] = pdb_file

    with open(yaml_file, "w") as file:
        yaml.dump(user_dict, file)

    main.run_platform_from_yaml(yaml_file)
