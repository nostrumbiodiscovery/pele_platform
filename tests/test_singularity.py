import os
import pytest
import pele_platform.constants.constants as cs
import pele_platform.Utilities.Helpers.yaml_parser as yp
import pele_platform.Checker.valid_flags as vf
from pele_platform.Utilities.Parameters.parameters import ParametersBuilder
from pele_platform.Utilities.Helpers.launcher import Launcher

test_path = os.path.join(cs.DIR, "Examples")

ARGS_1 = os.path.join(test_path, "singularity/input_adaptive.yaml")
ARGS_2 = os.path.join(test_path, "singularity/input_frag.yaml")

ADAPTIVE_PELE_EXEC = 'Pele_mpi'
ADAPTIVE_MPI_PARAMS = '"mpiParameters": "/path/to/singularity_container.sif",'
FRAG_PELE_EXEC = '/path/to/singularity_container.sif Pele_mpi'

testdata = [
    (ARGS_1, True, ADAPTIVE_PELE_EXEC, ADAPTIVE_MPI_PARAMS),
    (ARGS_2, False, FRAG_PELE_EXEC, None),
]

@pytest.mark.parametrize("ext_args, check_mpi, expected1, expected2", testdata)
def test_singularity_params(ext_args, check_mpi, expected1, expected2):
    """
    Test for the pele_exec and mpi_params values.

    Returns
    ----------
    boolean : result of the test.
    """
    # Parse yaml file
    yaml_obj = yp.YamlParser(ext_args, vf.VALID_FLAGS_PLATFORM)
    yaml_obj.read()

    # Prepare params
    launcher = Launcher(yaml_obj)
    launcher._define_package_to_run()
    params = launcher._args
    builder = ParametersBuilder()
    simulation_params = builder.build_adaptive_variables(params)

    errors = []

    if not simulation_params.pele_exec ==  expected1 :
        errors.append(
            f"conditions assert: [{simulation_params.pele_exec}] == [{expected1}] "
        )

    if check_mpi and not simulation_params.mpi_params ==  expected2:
        errors.append(
            f"conditions assert: [{simulation_params.mpi_params}] == [{expected2}] "
        )

    assert not errors

