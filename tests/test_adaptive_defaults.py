import os
import glob
import pele_platform.constants.constants as cs
import pele_platform.constants.pele_params as pp
import pele_platform.main as main

test_path = os.path.join(cs.DIR, "Examples")


BIAS_ARGS = os.path.join(test_path, "bias/input_defaults.yaml")
OUT_IN_ARGS = os.path.join(test_path, "out_in/input.yaml")
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

GLOBAL_DEFAULTS_ADAPTIVE = [
    '"type" : "inverselyProportional"',
    '"peleSteps" : 4,',
    '"iterations" : 100,',
    '"processors" : 250,',
    '[2.5, 5, 7]',
    '[1, 0.6, 0.0]'
]

GLOBAL_DEFAULTS_PELE = [
    '"numberOfStericTrials": 200',
    pp.GLOBAL
]

BIAS_DEFAULTS_ADAPTIVE = [
    '"type" : "epsilon"',
    '"peleSteps" : 8,',
    '"iterations" : 50,',
    '"processors" : 100,',
    '[1.5, 2, 5]',
    '[0.6, 0.4, 0.0]',
    '"metricColumnInReport" : 6',
    '"epsilon": 0.25'
]

BIAS_DEFAULTS_PELE = [
    '"numberOfStericTrials": 250',
    pp.BIAS
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
    pp.RESCORING
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

WATER_LIG_DEFAULTS_ADAPTIVE = [
    '[1.75, 2.5, 3.5, 5]',
    '"type" : "inverselyProportional"',
    '"peleSteps" : 12,',
    '"iterations" : 50',
    '[1.6, 1.2, 1, 0.0]'
]

WATER_DEFAULTS_ADAPTIVE = [
    '"type" : "independent"',
    '"peleSteps" : 12,',
    '"iterations" : 50',
    '"processors" : 128,'
]

WATER_LIG_DEFAULTS_PELE = [
    '"numberOfStericTrials": 100',
    '"overlapFactor": 0.5',
    '"radius": 10',
    pp.WATER_LIG
]

WATER_DEFAULTS_PELE = [
    pp.WATER_BS
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
    assert len(glob.glob(os.path.join(job.pele_dir, "input*.pdb"))) == (job.cpus-1)
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

def test_water_defaults(ext_args=WATER_ARGS):
    errors = []
    job = main.run_platform(ext_args)
    errors = check_file(job.pele_dir, "adaptive.conf", WATER_DEFAULTS_ADAPTIVE, errors)
    errors = check_file(job.pele_dir, "pele.conf", WATER_DEFAULTS_PELE, errors)
    assert not errors

def test_water_lig_defaults(ext_args=WATERLIG_ARGS):
    errors = []
    job = main.run_platform(ext_args)
    errors = check_file(job.pele_dir, "adaptive.conf", WATER_LIG_DEFAULTS_ADAPTIVE, errors)
    errors = check_file(job.pele_dir, "pele.conf", WATER_LIG_DEFAULTS_PELE, errors)
    assert not errors

def test_bias_defaults(ext_args=BIAS_ARGS):
    errors = []
    job = main.run_platform(ext_args)
    errors = check_file(job.pele_dir, "adaptive.conf", BIAS_DEFAULTS_ADAPTIVE, errors)
    errors = check_file(job.pele_dir, "pele.conf", BIAS_DEFAULTS_PELE, errors)
    assert not errors

def pca():
    pass

def out_in():
    pass

def test_rescoring_defaults(ext_args=RESCORING_ARGS):
    errors = []
    job = main.run_platform(ext_args)
    errors = check_file(job.pele_dir, "adaptive.conf", REF_DEFAULTS_ADAPTIVE, errors)
    errors = check_file(job.pele_dir, "pele.conf", REF_DEFAULTS_PELE, errors)
    assert not errors

def check_file(folder, filename, values, errors):
   filename = os.path.join(folder, filename)
   with open(filename, "r") as f:
      lines = f.readlines()
      for value in values:
          if value not in "".join(lines):
              errors.append(value) 
   return errors
