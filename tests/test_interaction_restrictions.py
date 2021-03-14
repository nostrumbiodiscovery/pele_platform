import os
import pytest
import pele_platform.constants.constants as cs
import pele_platform.Adaptive.interaction_restrictions as ir
import pele_platform.Utilities.Helpers.yaml_parser as yp
import pele_platform.Checker.valid_flags as vf
import pele_platform.main as main
from .test_adaptive import check_file

test_path = os.path.join(cs.DIR, "Examples")

ARGS_1 = os.path.join(test_path, "interaction_restrictions/input.yaml")
ARGS_2 = os.path.join(test_path, "interaction_restrictions/input_error.yaml")
ARGS_3 = os.path.join(test_path, "interaction_restrictions/input_error_2.yaml")
PDB_FILE = os.path.join(test_path, "interaction_restrictions/complex.pdb")

EXPECTED_METRICS = os.path.join(
    test_path, "interaction_restrictions/out/expected_metrics.txt"
)
EXPECTED_CONDITIONS = os.path.join(
    test_path, "interaction_restrictions/out/expected_conditions.txt"
)

INTERACTION_RESTRICTIONS_PELE = [
    '"interactionRestrictions":',
    '"distance0 < 3"',
    '"angle1 > 90"',
    '"angle1 < 180"',
]

METRIC_DISTANCE_PELE = [
    '"type":"com_distance"',
    '"tag":"distance0"',
    '"selection_group_1":{',
    '"atoms": { "ids":["A:318:_OG1"]}',
    '"selection_group_2":{',
    '"atoms": { "ids":["Z:201:_O3_"]}',
]

METRIC_ANGLE_PELE = [
    '"type":"atomsAngle"',
    '"tag":"angle1"',
    '"selection_group_1":{',
    '"atoms": { "ids":["A:318:_OG1"]}',
    '"selection_group_2":{',
    '"atoms": { "ids":["A:318:_HG1"]}',
    '"selection_group_3":{',
    '"atoms": { "ids":["Z:201:_O3_"]}',
]


def test_interaction_restrictions(ext_args=ARGS_1):
    """
    Test the pele.conf generated from a simulation with interaction restrictions.

    Parameters
    ----------
    ext_args : Path of the input.yaml file.

    Returns
    ----------
    boolean : result of the test.
    """
    errors = []
    job = main.run_platform(ext_args)
    errors = check_file(
        job.pele_dir, "pele.conf", INTERACTION_RESTRICTIONS_PELE, errors
    )
    errors = check_file(job.pele_dir, "pele.conf", METRIC_DISTANCE_PELE, errors)
    errors = check_file(job.pele_dir, "pele.conf", METRIC_ANGLE_PELE, errors)
    assert not errors


def test_metrics_and_conditions_to_json():
    """
    Test for the metrics_to_json and conditions_to_json methods.

    Returns
    ----------
    boolean : result of the test.
    """
    # Parse yaml file
    yaml_obj = yp.YamlParser(ARGS_1, vf.VALID_FLAGS_PLATFORM)
    yaml_obj.read()
    # Build interaction restrictions
    interaction_restrictions = ir.InteractionRestrictionsBuilder()
    interaction_restrictions.parse_interaction_restrictions(
        PDB_FILE, yaml_obj.interaction_restrictions
    )

    errors = []
    # Check metrics json output
    metrics = interaction_restrictions.metrics_to_json()
    expected_metrics_output = open(EXPECTED_METRICS, "r").read()
    if not metrics == expected_metrics_output:
        errors.append(f"Error in metrics json: {metrics} == {expected_metrics_output} ")

    # Check conditions json output
    conditions = interaction_restrictions.conditions_to_json()
    expected_conditions_output = open(EXPECTED_CONDITIONS, "r").read()
    if not conditions == expected_conditions_output:
        errors.append(
            f"conditions json assert: {conditions} == {expected_conditions_output} "
        )

    # assert no error message has been registered, else print messages
    assert not errors


"""
First test to throw an exception in a syntax error in input file (wrong number
of atoms).

Second test to throw an exception in a syntax error in input file (poorly defined metric).
"""


@pytest.mark.parametrize("file", [ARGS_2, ARGS_3])
def test_check_SyntaxError_exception(file):
    """
    Function to test an input file that should throw a Syntax Error exception.

    Parameters
    ----------
    file : Path of the input.yaml file.
    """
    # Parse yaml file
    yaml_obj = yp.YamlParser(file, vf.VALID_FLAGS_PLATFORM)
    yaml_obj.read()
    # Build interaction restrictions
    interaction_restrictions = ir.InteractionRestrictionsBuilder()
    with pytest.raises(SyntaxError):
        interaction_restrictions.parse_interaction_restrictions(
            PDB_FILE, yaml_obj.interaction_restrictions
        )
