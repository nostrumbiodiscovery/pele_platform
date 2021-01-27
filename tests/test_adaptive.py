import os
import pytest
import shutil
import pele_platform.constants.constants as cs
import pele_platform.main as main

test_path = os.path.join(cs.DIR, "Examples")

OUT_IN_ARGS = os.path.join(test_path, "out_in/input.yaml")
INDUCED_EX_ARGS = os.path.join(test_path, "induced_fit/input_exhaustive.yaml")
INDUCED_FAST_ARGS = os.path.join(test_path, "induced_fit/input_fast.yaml")
GLOBAL_ARGS = os.path.join(test_path, "global/input.yaml")
INPUTS_GLOBAL_ARGS = os.path.join(test_path, "global/input_inputs.yaml")
INPUTS_AST_ARGS = os.path.join(test_path, "global/input_inputs_asterisc.yaml")
EXIT_ARGS = os.path.join(test_path, "exit/input.yaml")
EXITSOFT_ARGS = os.path.join(test_path, "exit_soft/input.yaml")
WATER_ARGS = os.path.join(test_path, "water/input_bs.yaml")
RESTART_ARGS = os.path.join(test_path, "restart/input.yaml")
MSM_ARGS = os.path.join(test_path, "Msm/input.yaml")
MAE_ARGS = os.path.join(test_path, "induced_fit/input_mae.yaml")
PCA_ARGS = os.path.join(test_path, "pca/input.yaml")
PCA2_ARGS = os.path.join(test_path, "pca/input_str.yaml")
FLAGS_ARGS = os.path.join(test_path, "flags/input.yaml")
RESCORING_ARGS = os.path.join(test_path, "rescoring/input.yaml")
GPCR_ARGS = os.path.join(test_path, "gpcr/input.yaml")


ADAPTIVE_VALUES = [
    "water_processed.pdb",
    "SB4",
    '"outputPath": "output_sim"',
    '"processors" : 3',
    '"peleSteps" : 1,',
    '"iterations" : 1,',
    '"runEquilibration" : true,',
    '"equilibrationLength" : 11,',
    '"seed": 3000',
    '"values" : [1, 2, 3],',
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

PELE_VALUES = [
    "rep",
    "traj.xtc",
    "OBC",
    '"anmFrequency" : 3,',
    '"sideChainPredictionFrequency" : 3,',
    '"minimizationFrequency" : 3,',
    '"temperature": 3000,',
    '{ "type": "constrainAtomToPosition", "springConstant": 3, "equilibriumDistance": 0.0, "constrainThisAtom": "A:111:_CA_" },',
    '{ "type": "constrainAtomToPosition", "springConstant": 3, "equilibriumDistance": 0.0, "constrainThisAtom": "A:11:_CA_" }',
    '{ "type": "constrainAtomToPosition", "springConstant": 3, "equilibriumDistance": 0.0, "constrainThisAtom": "A:13:_CA_" }',
    '{ "type": "constrainAtomToPosition", "springConstant": 5, "equilibriumDistance": 0.0, "constrainThisAtom": "A:353:_CA_" }',
    '{ "type": "constrainAtomToPosition", "springConstant": 5, "equilibriumDistance": 0.0, "constrainThisAtom": "A:5:_CA_" }',
    '"radius": 3000',
    '"fixedCenter": [30,30,30]',
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
    "yaml",
    [
        INDUCED_EX_ARGS,
        INDUCED_FAST_ARGS,
        GLOBAL_ARGS,
        INPUTS_GLOBAL_ARGS,
        EXIT_ARGS,
        EXITSOFT_ARGS,
        OUT_IN_ARGS,
        RESTART_ARGS,
        RESCORING_ARGS,
        MAE_ARGS,
        INPUTS_AST_ARGS,
        GPCR_ARGS,
        PCA_ARGS,
    ],
)
def test_production_run(yaml):
    """
    Full PELE run in test mode (cpus = 5) using various YAML inputs to cover multiple simulation types.
    """
    main.run_platform(yaml)


def test_flags(ext_args=FLAGS_ARGS, output="NOR_solvent_OBC"):
    errors = []
    if os.path.exists(output):
        shutil.rmtree(output, ignore_errors=True)
    job = main.run_platform(ext_args)
    folder = job.folder
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
    if os.path.exists(os.path.join(folder, job.output)):
        errors.append("Debug flag not working")
    errors = check_file(folder, "adaptive.conf", ADAPTIVE_VALUES, errors)
    errors = check_file(folder, "pele.conf", PELE_VALUES, errors)
    errors = check_file(
        folder, "DataLocal/LigandRotamerLibs/SB4.rot.assign", "60", errors
    )
    assert not errors


@pytest.mark.parametrize("yaml", [PCA_ARGS, PCA2_ARGS])
def test_pca(yaml, output="PCA_result"):
    """
    Tests PCA runs in debug mode, both with an explicit list of PDB files for PCA and a wild card.
    """
    if os.path.exists(output):
        shutil.rmtree(output, ignore_errors=True)

    errors = []
    (job,) = main.run_platform(yaml)
    errors = check_file(job.pele_dir, "pele.conf", PCA_VALUES, errors)
    assert not errors


def check_file(folder, filename, values, errors):
    """
    Checks if file contains specific lines, otherwise adds missing values to errors.
    """
    filename = os.path.join(folder, filename)
    with open(filename, "r") as f:
        lines = f.readlines()
        for value in values:
            if value not in "".join(lines):
                errors.append(value)
    return errors
