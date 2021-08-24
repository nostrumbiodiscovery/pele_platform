import os
import pytest
import shutil
import yaml

import pele_platform.constants.constants as cs
import pele_platform.main as main
from tests.utils import check_file

test_path = os.path.join(cs.DIR, "Examples")


OUT_IN_ARGS = os.path.join(test_path, "out_in/input.yaml")
INDUCED_EX_ARGS = os.path.join(test_path, "induced_fit/input_exhaustive.yaml")
INDUCED_FAST_ARGS = os.path.join(test_path, "induced_fit/input_fast.yaml")
GLOBAL_ARGS = os.path.join(test_path, "global/input.yaml")
INPUTS_GLOBAL_ARGS = os.path.join(test_path, "global/input_inputs.yaml")
INPUTS_AST_ARGS = os.path.join(test_path, "global/input_inputs_asterisc.yaml")
EXIT_ARGS = os.path.join(test_path, "exit/input.yaml")
EXIT_SOFT_ARGS = os.path.join(test_path, "exit_soft/input.yaml")
MSM_ARGS = os.path.join(test_path, "Msm/input.yaml")
MAE_ARGS = os.path.join(test_path, "induced_fit/input_mae.yaml")
PCA_ARGS = os.path.join(test_path, "pca/input.yaml")
PCA2_ARGS = os.path.join(test_path, "pca/input_str.yaml")
FLAGS_ARGS = os.path.join(test_path, "flags/input.yaml")
RESCORING_ARGS = os.path.join(test_path, "rescoring/input.yaml")
GPCR_ARGS = os.path.join(test_path, "gpcr/input.yaml")
MINIMUM_STEPS_ARGS = os.path.join(test_path, "minimum_steps/input.yaml")

ADAPTIVE_VALUES = [
    "water_processed.pdb",
    "SB4",
    '"outputPath": "/home/agruzka/pele_platform/tests/NOR_solvent_OBC/output_sim",',
    '"processors" : 3',
    '"peleSteps" : 1,',
    '"iterations" : 1,',
    '"runEquilibration" : true,',
    '"equilibrationLength" : 11,',
    '"seed": 3000',
    '"values" : [1.0, 2.0, 3.0],',
    '"conditions": [0.1, 0.2, 0.3]',
    '"epsilon": 0.3',
    '"metricColumnInReport" : 3,',
    '"metricColumnInReport" : 3,',
    '"type" : "epsilon"',
    "exitContinuous",
    '"data": "done"',
    '"executable": "done"',
    '"documents": "done"',
]

GPCR_VALUES = [
    '"radius": 19.970223159033843,',
    '"fixedCenter": [-71.78435134887695,-13.431749963760375,-42.46209926605225]',
]

PELE_VALUES = [
    "rep",
    "traj.xtc",
    "OBC",
    '"anmFrequency" : 3,',
    '"sideChainPredictionFrequency" : 3,',
    '"minimizationFrequency" : 3,',
    '"temperature": 3000,',
    '{ "type": "constrainAtomToPosition", "springConstant": 3.0, "equilibriumDistance": 0.0, "constrainThisAtom": "A:111:_CA_" },',
    '{ "type": "constrainAtomToPosition", "springConstant": 3.0, "equilibriumDistance": 0.0, "constrainThisAtom": "A:11:_CA_" }',
    '{ "type": "constrainAtomToPosition", "springConstant": 3.0, "equilibriumDistance": 0.0, "constrainThisAtom": "A:13:_CA_" }',
    '{ "type": "constrainAtomToPosition", "springConstant": 5.0, "equilibriumDistance": 0.0, "constrainThisAtom": "A:353:_CA_" }',
    '{ "type": "constrainAtomToPosition", "springConstant": 5.0, "equilibriumDistance": 0.0, "constrainThisAtom": "A:5:_CA_" }',
    '"radius": 3000',
    '"fixedCenter": [30.0,30.0,30.0]',
    'tests/native.pdb"',
    '"atoms": { "ids":["A:6:_CG_"]}',
    '"atoms": { "ids":["A:6:_CD_"]}',
    "simulationLogPath",
    '"activateProximityDetection": false',
    '"overlapFactor": 3',
    '"steeringUpdateFrequency": 3',
    '"verboseMode": true,',
    '"displacementFactor" : 3',
    '"modesChangeFrequency" : 3,',
    '"directionGeneration" : "steered",',
    '"modesMixingOption" : "DontMixModes",',
    '"pickingCase" : "random"',
    '"numberOfModes": 3',
    '"preloadedModesIn" : "modes.nmd",',
    '"relaxationSpringConstant" : 3',
]

PCA_VALUES = [
    '"displacementFactor" : 2',
    '"modesChangeFrequency" : 50,',
    '"directionGeneration" : "oscillate",',
    '"modesMixingOption" : "doNotMixModes",',
    '"pickingCase" : "LOWEST_MODE"',
    '"numberOfModes": 1',
    '"relaxationSpringConstant" : 2',
]


@pytest.mark.parametrize(
    "yaml_file",
    [
        INDUCED_EX_ARGS,
        INDUCED_FAST_ARGS,
        GLOBAL_ARGS,
        INPUTS_GLOBAL_ARGS,
        EXIT_ARGS,
        EXIT_SOFT_ARGS,
        OUT_IN_ARGS,
        RESCORING_ARGS,
        MAE_ARGS,
        INPUTS_AST_ARGS,
        GPCR_ARGS,
        PCA_ARGS,
        PCA2_ARGS,
        MINIMUM_STEPS_ARGS,
    ],
)
def test_production_run(yaml_file):
    """
    Full PELE run in test mode (cpus = 5) using various YAML inputs to cover multiple simulation types.
    """
    main.run_platform_from_yaml(yaml_file)


@pytest.mark.parametrize("restart_type", ["restart", "adaptive_restart"])
def test_restart_flag(restart_type):
    """
    Checks if the platform can correctly restart and adaptive_restart a simulation from existing files (created in
    debug mode).
    """
    # First, run the platform in debug mode
    restart_yaml = os.path.join(test_path, "restart", "restart.yaml")
    job_parameters, = main.run_platform_from_yaml(restart_yaml)

    # Edit the original YAML file
    with open(restart_yaml, "r") as file:
        new_parameters = yaml.safe_load(file)

    new_parameters[restart_type] = True
    new_parameters["debug"] = False

    updated_restart_yaml = "updated_restart.yaml"
    with open(updated_restart_yaml, "w+") as new_file:
        yaml.dump(new_parameters, new_file)

    # Restart the simulation
    job_parameters2 = main.run_platform_from_yaml(updated_restart_yaml)
    print(job_parameters2)

    # Assert it reused existing directory and did not create a new one
    assert job_parameters.pele_dir == job_parameters2.pele_dir

    # Make sure pele.conf was not overwritten
    errors = check_file(
        job_parameters2.pele_dir, "pele.conf", 'overlapFactor": 0.5', []
    )
    assert not errors

    # Check if it finished
    assert os.path.exists(os.path.join(job_parameters2.pele_dir, "results"))


def test_flags(ext_args=FLAGS_ARGS, output="NOR_solvent_OBC"):
    errors = []
    if os.path.exists(output):
        shutil.rmtree(output, ignore_errors=True)
    builder, env = main.run_platform_from_yaml(ext_args)
    folder = env.folder
    if not os.path.exists(
        os.path.join(folder, "DataLocal/LigandRotamerLibs/STR.rot.assign")
    ) or not os.path.exists(
        os.path.join(folder, "DataLocal/LigandRotamerLibs/MG.rot.assign")
    ):
        errors.append("External rotamer flag not working")
    if not os.path.exists(
        os.path.join(folder, "DataLocal/Templates/OPLS2005/HeteroAtoms/strz")
    ) or not os.path.exists(
        os.path.join(folder, "DataLocal/Templates/OPLS2005/HeteroAtoms/mgz")
    ):
        errors.append("External templates flag not working")
    if os.path.exists(os.path.join(folder, env.output)):
        errors.append("Debug flag not working")
    errors = check_file(folder, "adaptive.conf", ADAPTIVE_VALUES, errors)
    errors = check_file(folder, "pele.conf", PELE_VALUES, errors)
    errors = check_file(
        folder, "DataLocal/LigandRotamerLibs/SB4.rot.assign", "60", errors
    )
    assert not errors
