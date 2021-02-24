import os
import glob
import pele_platform.constants.constants as cs
import pele_platform.constants.pele_params as pp
import pele_platform.main as main

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
    '"processors" : 60,'
    ]

INDUCE_FIT_FAST_DEFAULTS_ADAPTIVE = [
    '"type" : "inverselyProportional"',
    '"iterations" : 30,',
    '"peleSteps" : 12,',
    '"processors" : 60,'
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
    pp.OUT_IN
]

REF_DEFAULTS_ADAPTIVE = [
    '"type" : "independent"',
    '"iterations" : 20,',
    '"peleSteps" : 12,',
    '"processors" : 60,'
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
    '"type": "exitContinuous"'
]

EXIT_SOFT_DEFAULTS_ADAPTIVE = [
    '"type" : "independentMetric"',
    '"iterations" : 1000',
    '"processors" : 25',
    '[1, 2.5]',
    '[1.1]',
    '"peleSteps" : 2,',
    '"type": "exitContinuous"',
    '"condition": "max"'
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
     '"fixedCenter": [-71.78435134887695,-13.431749963760375,-42.46209926605225]',
     '"numberOfStericTrials": 100,',
     pp.GPCR_ORTH
]

GPCR_DEFAULTS_ADAPTIVE = [
    '"type" : "epsilon"',
    '"metricColumnInReport" : 6',
    '"epsilon": 0.25',
    '"iterations" : 50',
    '"peleSteps" : 8',
    '"values" : [1.75, 2.5, 4],',
    '"conditions": [0.7, 0.4, 0.0]'
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


def test_induced_exhaustive_defaults(ext_args=INDUCED_EX_ARGS):
    errors = []
    job = main.run_platform(ext_args)
    errors = check_file(job.pele_dir, "adaptive.conf", INDUCE_FIT_EXHAUSTIVE_DEFAULTS_ADAPTIVE, errors)
    errors = check_file(job.pele_dir, "pele.conf", INDUCE_FIT_PELE, errors)
    assert not errors

def test_induced_fast_defaults(ext_args=INDUCED_FAST_ARGS):
    errors = []
    job = main.run_platform(ext_args)
    errors = check_file(job.pele_dir, "adaptive.conf", INDUCE_FIT_FAST_DEFAULTS_ADAPTIVE, errors)
    errors = check_file(job.pele_dir, "pele.conf", INDUCE_FIT_PELE, errors)
    assert not errors

def test_global_defaults(ext_args=GLOBAL_ARGS):
    errors = []
    job = main.run_platform(ext_args)
    errors = check_file(job.pele_dir, "adaptive.conf", GLOBAL_DEFAULTS_ADAPTIVE, errors)
    errors = check_file(job.pele_dir, "pele.conf", GLOBAL_DEFAULTS_PELE, errors)
    assert len(glob.glob(os.path.join(job.inputs_dir, "input*.pdb"))) == (job.cpus-1)
    assert not errors

def test_exit_defaults(ext_args=EXIT_ARGS):
    errors = []
    job = main.run_platform(ext_args)
    errors = check_file(job.pele_dir, "adaptive.conf", EXIT_DEFAULTS_ADAPTIVE, errors)
    errors = check_file(job.pele_dir, "pele.conf", EXIT_DEFAULTS_PELE, errors)
    assert not errors

def test_softexit_defaults(ext_args=EXITSOFT_ARGS):
    errors = []
    job = main.run_platform(ext_args)
    errors = check_file(job.pele_dir, "adaptive.conf", EXIT_SOFT_DEFAULTS_ADAPTIVE, errors)
    errors = check_file(job.pele_dir, "pele.conf", EXIT_DEFAULTS_PELE, errors)
    assert not errors

def test_water_lig_defaults(ext_args=WATERLIG_ARGS):
    errors = []
    job = main.run_platform(ext_args)
    errors = check_file(job.pele_dir, "pele.conf", WATER_PARAMS_DEFAULTS_PELE, errors)
    assert not errors

def test_out_in_defaults(ext_args=OUT_IN_ARGS):
    errors = []
    job = main.run_platform(ext_args)
    errors = check_file(job.pele_dir, "adaptive.conf", OUT_IN_DEFAULTS_ADAPTIVE, errors)
    errors = check_file(job.pele_dir, "pele.conf", OUT_IN_DEFAULTS_PELE, errors)
    assert not errors

def test_bias_defaults(ext_args=BIAS_ARGS):
    errors = []
    job = main.run_platform(ext_args)
    errors = check_file(job.pele_dir, "adaptive.conf", BIAS_DEFAULTS_ADAPTIVE, errors)
    assert not errors

def test_rescoring_defaults(ext_args=RESCORING_ARGS):
    errors = []
    job = main.run_platform(ext_args)
    errors = check_file(job.pele_dir, "adaptive.conf", REF_DEFAULTS_ADAPTIVE, errors)
    errors = check_file(job.pele_dir, "pele.conf", REF_DEFAULTS_PELE, errors)
    assert not errors

def test_gpcr_defaults(ext_args=GPCR_ARGS):
    errors = []
    job = main.run_platform(ext_args)
    errors = check_file(job.pele_dir, "adaptive.conf", GPCR_DEFAULTS_ADAPTIVE, errors)
    errors = check_file(job.pele_dir, "pele.conf", GPCR_DEFAULTS_PELE, errors)
    assert not errors

def test_gpcr2_defaults(ext_args=GPCR2_ARGS):
    errors = []
    job = main.run_platform(ext_args)
    errors = check_file(job.pele_dir, "adaptive.conf", GPCR_DEFAULTS_ADAPTIVE, errors)
    errors = check_file(job.pele_dir, "pele.conf", GPCR2_DEFAULTS_PELE, errors)
    assert not errors

def check_file(folder, filename, values, errors):
   filename = os.path.join(folder, filename)
   with open(filename, "r") as f:
      lines = f.readlines()
      for value in values:
          if value not in "".join(lines):
              errors.append(value) 
   return errors
