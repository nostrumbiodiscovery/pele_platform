import pytest
import os
import pele_platform.constants.constants as cs
import pele_platform.main as main

test_path = os.path.join(cs.DIR, "Examples")


BIAS_ARGS = [os.path.join(test_path, "bias/input.yaml")]
OUT_IN_ARGS = [os.path.join(test_path, "out_in/input.yaml")]
INDUCED_ARGS = [os.path.join(test_path, "induced_fit/input.yaml")]
GLOBAL_ARGS = [os.path.join(test_path, "global/input.yaml")]
EXIT_ARGS = [os.path.join(test_path, "exit/input.yaml")]
EXITSOFT_ARGS = [os.path.join(test_path, "exit_soft/input.yaml")]
WATER_ARGS = [os.path.join(test_path, "water/input_bs.yaml")]
WATERLIG_ARGS = [os.path.join(test_path, "water/input_lig.yaml")]
MSM_ARGS = [os.path.join(test_path, "Msm/input.yaml")]



def test_out_in(ext_args=OUT_IN_ARGS):
    arguments = main.parseargs_yaml(ext_args)
    arguments = main.YamlParser(arguments.input_file)
    main.set_software_to_use(arguments)
    main.Launcher(arguments).launch()

def test_induced(ext_args=INDUCED_ARGS):
    arguments = main.parseargs_yaml(ext_args)
    arguments = main.YamlParser(arguments.input_file)
    main.set_software_to_use(arguments)
    main.Launcher(arguments).launch()

def test_exit(ext_args=EXIT_ARGS):
    arguments = main.parseargs_yaml(ext_args)
    arguments = main.YamlParser(arguments.input_file)
    main.set_software_to_use(arguments)
    main.Launcher(arguments).launch()

def test_exitsoft(ext_args=EXITSOFT_ARGS):
    arguments = main.parseargs_yaml(ext_args)
    arguments = main.YamlParser(arguments.input_file)
    main.set_software_to_use(arguments)
    main.Launcher(arguments).launch()

def test_water(ext_args=WATER_ARGS):
    arguments = main.parseargs_yaml(ext_args)
    arguments = main.YamlParser(arguments.input_file)
    main.set_software_to_use(arguments)
    main.Launcher(arguments).launch()

def test_water_lig(ext_args=WATERLIG_ARGS):
    arguments = main.parseargs_yaml(ext_args)
    arguments = main.YamlParser(arguments.input_file)
    main.set_software_to_use(arguments)
    main.Launcher(arguments).launch()

def test_bias(ext_args=BIAS_ARGS):
    arguments = main.parseargs_yaml(ext_args)
    arguments = main.YamlParser(arguments.input_file)
    main.set_software_to_use(arguments)
    main.Launcher(arguments).launch()

def test_msm(ext_args=MSM_ARGS):
    arguments = main.parseargs_yaml(ext_args)
    arguments = main.YamlParser(arguments.input_file)
    main.set_software_to_use(arguments)
    main.Launcher(arguments).launch()

