import os
import pytest

import pele_platform.constants.constants as constants
import pele_platform.Utilities.Helpers.constraints as constraints
from pele_platform import main
from . import test_adaptive as ta

test_path = os.path.join(constants.DIR, "Examples")

terminal_lines_raw = [
    '{ "type": "constrainAtomToPosition", "springConstant": 5, "equilibriumDistance": 0.0, "constrainThisAtom": "A:1:_CA_" },',
    '{ "type": "constrainAtomToPosition", "springConstant": 5, "equilibriumDistance": 0.0, "constrainThisAtom": "A:322:_CA_" },',
    '{ "type": "constrainAtomToPosition", "springConstant": 5, "equilibriumDistance": 0.0, "constrainThisAtom": "B:1:_CA_" },',
    '{ "type": "constrainAtomToPosition", "springConstant": 5, "equilibriumDistance": 0.0, "constrainThisAtom": "B:171:_CA_" },',
]

terminal_lines_file = [
    '{ "type": "constrainAtomToPosition", "springConstant": 5, "equilibriumDistance": 0.0, "constrainThisAtom": "A:1:_CA_" },',
    '{ "type": "constrainAtomToPosition", "springConstant": 5, "equilibriumDistance": 0.0, "constrainThisAtom": "A:322:_CA_" },',
    '{ "type": "constrainAtomToPosition", "springConstant": 5, "equilibriumDistance": 0.0, "constrainThisAtom": "B:1:_CA_" },',
    '{ "type": "constrainAtomToPosition", "springConstant": 5, "equilibriumDistance": 0.0, "constrainThisAtom": "B:171:_CA_" }',
]

default_interval_lines = [
    '{ "type": "constrainAtomToPosition", "springConstant": 0.5, "equilibriumDistance": 0.0, "constrainThisAtom": "A:11:_CA_" },',
    '{ "type": "constrainAtomToPosition", "springConstant": 0.5, "equilibriumDistance": 0.0, "constrainThisAtom": "A:21:_CA_" },',
    '{ "type": "constrainAtomToPosition", "springConstant": 0.5, "equilibriumDistance": 0.0, "constrainThisAtom": "B:11:_CA_" },',
    '{ "type": "constrainAtomToPosition", "springConstant": 0.5, "equilibriumDistance": 0.0, "constrainThisAtom": "B:21:_CA_" },',
]

interval7_lines = [
    '{ "type": "constrainAtomToPosition", "springConstant": 0.7, "equilibriumDistance": 0.0, "constrainThisAtom": "A:8:_CA_" },',
    '{ "type": "constrainAtomToPosition", "springConstant": 0.7, "equilibriumDistance": 0.0, "constrainThisAtom": "A:15:_CA_" },',
    '{ "type": "constrainAtomToPosition", "springConstant": 0.7, "equilibriumDistance": 0.0, "constrainThisAtom": "B:8:_CA_" },',
    '{ "type": "constrainAtomToPosition", "springConstant": 0.7, "equilibriumDistance": 0.0, "constrainThisAtom": "B:15:_CA_" },',
]

custom_lines = [
    '{ "type": "constrainAtomToPosition", "springConstant": 0.7, "equilibriumDistance": 0.0, "constrainThisAtom": "A:16:_CA_" },',
    '{ "type": "constrainAtomToPosition", "springConstant": 0.7, "equilibriumDistance": 0.0, "constrainThisAtom": "A:31:_CA_" },',
    '{ "type": "constrainAtomToPosition", "springConstant": 0.7, "equilibriumDistance": 0.0, "constrainThisAtom": "B:16:_CA_" },',
    '{ "type": "constrainAtomToPosition", "springConstant": 0.7, "equilibriumDistance": 0.0, "constrainThisAtom": "B:31:_CA_" },',
]

custom_ter_lines = [
    '{ "type": "constrainAtomToPosition", "springConstant": 77, "equilibriumDistance": 0.0, "constrainThisAtom": "A:1:_CA_" },',
    '{ "type": "constrainAtomToPosition", "springConstant": 77, "equilibriumDistance": 0.0, "constrainThisAtom": "A:322:_CA_" },',
    '{ "type": "constrainAtomToPosition", "springConstant": 77, "equilibriumDistance": 0.0, "constrainThisAtom": "B:1:_CA_" },',
    '{ "type": "constrainAtomToPosition", "springConstant": 77, "equilibriumDistance": 0.0, "constrainThisAtom": "B:171:_CA_" },',
]


@pytest.mark.parametrize(
    "input_yaml",
    [
        os.path.join(test_path, "constraints/input_ca_noppp.yaml"),
        os.path.join(test_path, "constraints/input_ca_ppp.yaml"),
    ],
)
def test_ca_constraints_production(input_yaml):
    """
    Tests backbone and terminal constraints in debug mode, both with and without PPP preprocessing.
    """
    errors = []
    job = main.run_platform(input_yaml)
    errors = ta.check_file(
        job.pele_dir, "pele.conf", terminal_lines_file + default_interval_lines, errors
    )
    print("errors\n", errors)
    assert not errors


@pytest.mark.parametrize(
    (
        "interval",
        "backbone_spring",
        "terminal_spring",
        "expected_ca_lines",
        "expected_ter_lines",
    ),
    [
        (7, 0.7, 5, interval7_lines, terminal_lines_raw),
        (15, 0.7, 77, custom_lines, custom_ter_lines),
    ],
)
def test_ca_constraints_builder(
    interval, backbone_spring, terminal_spring, expected_ca_lines, expected_ter_lines
):
    """
    Test the alpha carbon constraint builder with various parameters.
    """
    obj = constraints.ConstraintBuilder(
        os.path.join(test_path, "constraints/ca_constraints.pdb"),
        interval=interval,
        backbone_spring=backbone_spring,
        terminal_spring=terminal_spring,
    )
    obj.create_constraints()

    assert obj.terminal_constraints == expected_ter_lines

    for line in expected_ca_lines:
        assert line in obj.backbone_constraints


def test_find_gaps():
    # TODO: Find a system with gaps and test the prody/PPP function.
    pass
