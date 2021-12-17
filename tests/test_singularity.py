import os
import pytest
import pele_platform.constants.constants as cs
import pele_platform.Utilities.Helpers.yaml_parser as yp
import pele_platform.Checker.valid_flags as vf
from pele_platform.Utilities.Parameters.parameters import ParametersBuilder
from pele_platform.Utilities.Helpers.launcher import Launcher
import pele_platform.Checker.main as ck
import pele_platform.Errors.custom_errors as ce
from pele_platform import main

test_path = os.path.join(cs.DIR, "Examples")

ARGS_1 = os.path.join(test_path, "singularity/input_adaptive.yaml")
ARGS_2 = os.path.join(test_path, "singularity/input_frag.yaml")

ADAPTIVE_PELE_EXEC = "Pele_mpi"
ADAPTIVE_MPI_PARAMS = '"/path/to/singularity_container.sif",'
FRAG_PELE_EXEC = "/path/to/singularity_container.sif Pele_mpi"

testdata = [
    (ARGS_1, True, ADAPTIVE_PELE_EXEC, ADAPTIVE_MPI_PARAMS),
    (ARGS_2, False, FRAG_PELE_EXEC, None),
]


def test_singularity_exec(ext_args=ARGS_1):
    """
    Test singularity checker (throw an ExecutableNotInPath exception).

    Parameters
    ----------
    ext_args : Path of the input.yaml file.

    Returns
    ----------
    boolean : result of the test.
    """
    # Prepare params
    params = _read_args(ext_args)

    with pytest.raises(ce.ExecutableNotInPath):
        ck.check_executable_and_env_variables(params)


@pytest.mark.parametrize("ext_args, check_mpi, expected1, expected2", testdata)
def test_singularity_params(ext_args, check_mpi, expected1, expected2):
    """
    Test for the pele_exec and mpi_params values.

    Parameters
    ----------
    ext_args : Path of the input.yaml file.
    check_mpi: Boolean, activates mpi_params check.
    expected1: Expected value for pele_exec.
    expected2: Expected value for mpi_params.

    Returns
    ----------
    boolean : result of the test.
    """
    # Prepare params
    params = _read_args(ext_args)
    builder = ParametersBuilder()
    simulation_params = builder.build_adaptive_variables(params)

    # Check parameters
    assert simulation_params.pele_exec == expected1
    if check_mpi:
        mpi_params_name = "srunParameters" if params.usesrun else "mpiParameters"
        mpi_expected_params = f'"{mpi_params_name}": {expected2}'
        assert simulation_params.pele_mpi_params == mpi_expected_params


def test_singularity_mpi_params():
    """
    Checks if custom pele_mpi_parameters are correctly set when using singularity.
    """
    yaml_file = os.path.join(test_path, "singularity/input_adaptive_mpi_params.yaml")
    job = main.run_platform_from_yaml(yaml_file)

    if job.usesrun == "false":
        params_name = "mpiParameters"
    else:
        params_name = "srunParameters"
    assert job.pele_mpi_params == f'"{params_name}": "--prefix custom parameter /path/to/singularity_container.sif",'


def _read_args(file):
    """
    Internal function: parse yaml file and prepare args.

    Parameters
    ----------
    file : str
        Path of the input.yaml file.

    Returns
    ----------
    args : simulation args.
    """
    # Parse yaml file
    yaml_obj = yp.YamlParser(file, vf.VALID_FLAGS_PLATFORM)
    yaml_obj.read()

    # Prepare params
    launcher = Launcher(yaml_obj)
    launcher._define_package_to_run()

    return launcher._args
