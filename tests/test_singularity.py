import os
import pytest
import pele_platform.constants.constants as cs
import pele_platform.Utilities.Helpers.yaml_parser as yp
from pele_platform.Utilities.Parameters.parameters import ParametersBuilder
from pele_platform.Utilities.Helpers.launcher import Launcher
import pele_platform.Errors.custom_errors as ce

test_path = os.path.join(cs.DIR, "Examples")

ARGS_1 = os.path.join(test_path, "singularity/input_adaptive.yaml")
ARGS_2 = os.path.join(test_path, "singularity/input_frag.yaml")

ADAPTIVE_PELE_EXEC = "Pele_mpi"
ADAPTIVE_MPI_PARAMS = '"/path/to/singularity_container.sif",'
FRAG_PELE_EXEC = "/path/to/singularity_container.sif Pele_mpi"

test_data = [
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
    yaml_obj = yp.YamlParser(ext_args)
    yaml_obj.read()
    launcher = Launcher(yaml_obj)

    with pytest.raises(ce.ExecutableNotInPath):
        launcher.launch()


@pytest.mark.parametrize("ext_args, check_mpi, expected1, expected2", test_data)
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
        assert simulation_params.mpi_params == mpi_expected_params


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
    yaml_obj = yp.YamlParser(file)
    yaml_obj.read()

    # Prepare params
    launcher = Launcher(yaml_obj)
    launcher.launch()

    return launcher._args
