import os
import pytest
import pele_platform.constants.constants as cs
import pele_platform.constants.pele_params as pp
import pele_platform.main as main
from . import test_adaptive as tk

test_path = os.path.join(cs.DIR, "Examples")

OUT_IN_ARGS = os.path.join(test_path, "out_in/input_default.yaml")
BIAS_ARGS = os.path.join(test_path, "bias/input_defaults.yaml")
INDUCED_EX_ARGS = os.path.join(test_path, "induced_fit/input_exhaustive_defaults.yaml")
INDUCED_FAST_ARGS = os.path.join(test_path, "induced_fit/input_fast_defaults.yaml")
GLOBAL_ARGS = os.path.join(test_path, "global/input_defaults.yaml")
EXIT_ARGS = os.path.join(test_path, "exit/input_defaults.yaml")
EXITSOFT_ARGS = os.path.join(test_path, "exit_soft/input_defaults.yaml")
WATER_ARGS = os.path.join(test_path, "water/input_bs_defaults.yaml")
ALL_WATER_ARGS = os.path.join(test_path, "water/input_all.yaml")
WATERLIG_ARGS = os.path.join(test_path, "water/input_lig_defaults.yaml")
RESTART_ARGS = os.path.join(test_path, "restart/input.yaml")
MSM_ARGS = os.path.join(test_path, "Msm/input.yaml")
MAE_ARGS = os.path.join(test_path, "induced_fit/input_mae.yaml")
PCA_ARGS = os.path.join(test_path, "pca/input.yaml")
FLAGS_ARGS = os.path.join(test_path, "flags/input.yaml")
RESCORING_ARGS = os.path.join(test_path, "rescoring/input_defaults.yaml")
GPCR_ARGS = os.path.join(test_path, "gpcr/input_defaults.yaml")
GPCR2_ARGS = os.path.join(test_path, "gpcr/input_defaults2.yaml")

INDUCE_FIT_EXHAUSTIVE_DEFAULTS_ADAPTIVE = [
    '"type" : "independent"',
    '"iterations" : 1,',
    '"peleSteps" : 1000,',
    '"processors" : 60,',
]

INDUCE_FIT_FAST_DEFAULTS_ADAPTIVE = [
    '"type" : "inverselyProportional"',
    '"iterations" : 30,',
    '"peleSteps" : 12,',
    '"processors" : 60,',
]

INDUCE_FIT_PELE = [
    pp.INDUCED_FIT,
    '"radius": 6,',
    '"numberOfStericTrials": 500'
    ]

BIAS_DEFAULTS_ADAPTIVE = [
    '"type" : "epsilon"',
    '"epsilon": 1'
]

GLOBAL_DEFAULTS_ADAPTIVE = [
    '"type" : "inverselyProportional"',
    '"peleSteps" : 8,',
    '"iterations" : 100,',
    '"processors" : 250,',
    '[2.0, 5, 7]',
    '[1, 0.6, 0.0]'
]

GLOBAL_DEFAULTS_PELE = [
    '"numberOfStericTrials": 200',
    pp.GLOBAL
]

OUT_IN_DEFAULTS_ADAPTIVE = [
    '"type" : "inverselyProportional"',
    '"peleSteps" : 8,',
    '"iterations" : 100,',
    '[2, 5, 7]',
    '[1, 0.6, 0.0]',
]

OUT_IN_DEFAULTS_PELE = [
    '"numberOfStericTrials": 250',
    '"overlapFactor": 0.65',
    pp.OUT_IN,
]

REF_DEFAULTS_ADAPTIVE = [
    '"type" : "independent"',
    '"iterations" : 20,',
    '"peleSteps" : 12,',
    '"processors" : 60,',
]

REF_DEFAULTS_PELE = [
    '"numberOfStericTrials": 500',
    '"radius": 6,',
    '"anmFrequency" : 6,',
    '"sideChainPredictionFrequency" : 3,',
    '"minimizationFrequency" : 1,',
    '"temperature": 1000',
    '"displacementFactor" : 0.5',
    '"modesChangeFrequency" : 3,',
    pp.RESCORING,
    pp.RESCORING,
    '{ "type": "constrainAtomToPosition", "springConstant": 5.0, "equilibriumDistance": 0.0, "constrainThisAtom": "A:1:_CA_" },',
    '{ "type": "constrainAtomToPosition", "springConstant": 5.0, "equilibriumDistance": 0.0, "constrainThisAtom": "A:251:_CA_" }',
    '{ "type": "constrainAtomToPosition", "springConstant": 2.5, "equilibriumDistance": 0.0, "constrainThisAtom": "A:9:_CA_" }',
    '{ "type": "constrainAtomToPosition", "springConstant": 2.5, "equilibriumDistance": 0.0, "constrainThisAtom": "A:17:_CA_" }',
]

EXIT_DEFAULTS_ADAPTIVE = [
    '"type" : "epsilon"',
    '"metricColumnInReport" : 6',
    '"epsilon": 0.75',
    '"iterations" : 1000',
    '"processors" : 25',
    '[1, 2.5]',
    '[1.1]',
    '"peleSteps" : 2,',
    '"type": "exitContinuous"',
]

EXIT_SOFT_DEFAULTS_ADAPTIVE = [
    '"type" : "independentMetric"',
    '"iterations" : 1000',
    '"processors" : 25',
    '[1, 2.5]',
    '[1.1]',
    '"peleSteps" : 2,',
    '"type": "exitContinuous"',
    '"condition": "max"',
]


EXIT_DEFAULTS_PELE = [
    '"radius": 10',
    '"numberOfStericTrials": 500',
    pp.IN_OUT

]

WATER_PARAMS_DEFAULTS_PELE = [
    pp.WATER_PARAMS,
    '"watersToPerturb": {"links": {"ids": ["M:1"] }}' 
]

GPCR_DEFAULTS_PELE = [
    '"radius": 19.970223159033843,',
    '"fixedCenter": [-71.78435134887695, -13.431749963760375, -42.46209926605225]',
    '"numberOfStericTrials": 100,',
    pp.GPCR_ORTH,
]

GPCR_DEFAULTS_ADAPTIVE = [
    '"type" : "epsilon"',
    '"metricColumnInReport" : 6',
    '"epsilon": 0.25',
    '"iterations" : 50',
    '"peleSteps" : 8',
    '"values" : [1.75, 2.5, 4],',
    '"conditions": [0.7, 0.4, 0.0]',
]

GPCR2_DEFAULTS_PELE = [
     '"radius": 10,',
     '"fixedCenter": [10,10,10]',
     '"numberOfStericTrials": 100,',
     pp.GPCR_ORTH,
    '{ "type": "constrainAtomToPosition", "springConstant": 5.0, "equilibriumDistance": 0.0, "constrainThisAtom": "A:65:_CA_" },',
    '{ "type": "constrainAtomToPosition", "springConstant": 5.0, "equilibriumDistance": 0.0, "constrainThisAtom": "A:70:_CA_" },',
    '{ "type": "constrainAtomToPosition", "springConstant": 5.0, "equilibriumDistance": 0.0, "constrainThisAtom": "A:60:_CA_" },',
    '{ "type": "constrainAtomToPosition", "springConstant": 5.0, "equilibriumDistance": 0.0, "constrainThisAtom": "A:347:_CA_" }'
]


@pytest.mark.parametrize(
    ("yaml", "expected_adaptive", "expected_pele"),
    [
        (INDUCED_EX_ARGS, INDUCE_FIT_EXHAUSTIVE_DEFAULTS_ADAPTIVE, INDUCE_FIT_PELE),
        (INDUCED_FAST_ARGS, INDUCE_FIT_FAST_DEFAULTS_ADAPTIVE, INDUCE_FIT_PELE),
        (GLOBAL_ARGS, GLOBAL_DEFAULTS_ADAPTIVE, GLOBAL_DEFAULTS_PELE),
        (EXIT_ARGS, EXIT_DEFAULTS_ADAPTIVE, EXIT_DEFAULTS_PELE),
        (EXITSOFT_ARGS, EXIT_SOFT_DEFAULTS_ADAPTIVE, EXIT_DEFAULTS_PELE),
        (OUT_IN_ARGS, OUT_IN_DEFAULTS_ADAPTIVE, OUT_IN_DEFAULTS_PELE),
        (WATERLIG_ARGS, "", WATER_PARAMS_DEFAULTS_PELE),
        (BIAS_ARGS, BIAS_DEFAULTS_ADAPTIVE, ""),
        (RESCORING_ARGS, REF_DEFAULTS_ADAPTIVE, REF_DEFAULTS_PELE),
        (GPCR_ARGS, GPCR_DEFAULTS_ADAPTIVE, GPCR_DEFAULTS_PELE),
        (GPCR2_ARGS, GPCR_DEFAULTS_ADAPTIVE, GPCR2_DEFAULTS_PELE),
    ],
)
def test_all_defaults(yaml, expected_adaptive, expected_pele):
    """
    Runs all simulation types in debug mode to ensure the correct default values are set in adaptive.conf and pele.conf.
    """
    errors = []
    output = main.run_platform_from_yaml(yaml)
    folder = output[0].pele_dir if type(output) == list else output.pele_dir
    errors = tk.check_file(folder, "adaptive.conf", expected_adaptive, errors)
    errors = tk.check_file(folder, "pele.conf", expected_pele, errors)
    assert not errors
