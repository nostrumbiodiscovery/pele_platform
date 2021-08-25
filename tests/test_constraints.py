import os
import pytest
import yaml

import pele_platform.constants.constants as constants
import pele_platform.Utilities.Helpers.constraints.alpha_constraints as alpha_constraints
import tests.utils
from pele_platform import main
from pele_platform.Utilities.Helpers.constraints import smiles_constraints as smi
from tests.test_others import EXTERNAL_CONSTR_ARGS, EXT_CONSTR, ARGS_SMARTS_CONSTR

test_path = os.path.join(constants.DIR, "Examples")

terminal_lines_raw = [
    '{ "type": "constrainAtomToPosition", "springConstant": 5.0, "equilibriumDistance": 0.0, "constrainThisAtom": "A:1:_CA_" },',
    '{ "type": "constrainAtomToPosition", "springConstant": 5.0, "equilibriumDistance": 0.0, "constrainThisAtom": "A:322:_CA_" },',
    '{ "type": "constrainAtomToPosition", "springConstant": 5.0, "equilibriumDistance": 0.0, "constrainThisAtom": "B:1:_CA_" },',
    '{ "type": "constrainAtomToPosition", "springConstant": 5.0, "equilibriumDistance": 0.0, "constrainThisAtom": "B:171:_CA_" },',
]

terminal_lines_file = [
    '{ "type": "constrainAtomToPosition", "springConstant": 5.0, "equilibriumDistance": 0.0, "constrainThisAtom": "A:1:_CA_" },',
    '{ "type": "constrainAtomToPosition", "springConstant": 5.0, "equilibriumDistance": 0.0, "constrainThisAtom": "A:322:_CA_" },',
    '{ "type": "constrainAtomToPosition", "springConstant": 5.0, "equilibriumDistance": 0.0, "constrainThisAtom": "B:1:_CA_" },',
    '{ "type": "constrainAtomToPosition", "springConstant": 5.0, "equilibriumDistance": 0.0, "constrainThisAtom": "B:171:_CA_" }',
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
    '{ "type": "constrainAtomToPosition", "springConstant": 77.0, "equilibriumDistance": 0.0, "constrainThisAtom": "A:1:_CA_" },',
    '{ "type": "constrainAtomToPosition", "springConstant": 77.0, "equilibriumDistance": 0.0, "constrainThisAtom": "A:322:_CA_" },',
    '{ "type": "constrainAtomToPosition", "springConstant": 77.0, "equilibriumDistance": 0.0, "constrainThisAtom": "B:1:_CA_" },',
    '{ "type": "constrainAtomToPosition", "springConstant": 77.0, "equilibriumDistance": 0.0, "constrainThisAtom": "B:171:_CA_" },',
]

set1_lines = [
    '{ "type": "constrainAtomToPosition", "springConstant": 2.5, "equilibriumDistance": 0.0, "constrainThisAtom": "A:8:_CA_" },',
    '{ "type": "constrainAtomToPosition", "springConstant": 2.5, "equilibriumDistance": 0.0, "constrainThisAtom": "A:15:_CA_" },',
    '{ "type": "constrainAtomToPosition", "springConstant": 5.0, "equilibriumDistance": 0.0, "constrainThisAtom": "A:1:_CA_" },',
    '{ "type": "constrainAtomToPosition", "springConstant": 5.0, "equilibriumDistance": 0.0, "constrainThisAtom": "A:251:_CA_" }',
]

set2_lines = [
    '{ "type": "constrainAtomToPosition", "springConstant": 0.7, "equilibriumDistance": 0.0, "constrainThisAtom": "A:79:_CA_" },',
    '{ "type": "constrainAtomToPosition", "springConstant": 0.7, "equilibriumDistance": 0.0, "constrainThisAtom": "A:156:_CA_" },',
    '{ "type": "constrainAtomToPosition", "springConstant": 7.0, "equilibriumDistance": 0.0, "constrainThisAtom": "A:1:_CA_" },',
    '{ "type": "constrainAtomToPosition", "springConstant": 7.0, "equilibriumDistance": 0.0, "constrainThisAtom": "A:251:_CA_" }',
]

set3_lines = [
    '{ "type": "constrainAtomToPosition", "springConstant": 11.0, "equilibriumDistance": 0.0, "constrainThisAtom": "A:12:_CA_" },',
    '{ "type": "constrainAtomToPosition", "springConstant": 11.0, "equilibriumDistance": 0.0, "constrainThisAtom": "A:23:_CA_" },',
    '{ "type": "constrainAtomToPosition", "springConstant": 5.0, "equilibriumDistance": 0.0, "constrainThisAtom": "A:1:_CA_" },',
    '{ "type": "constrainAtomToPosition", "springConstant": 5.0, "equilibriumDistance": 0.0, "constrainThisAtom": "A:251:_CA_" }',
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
    Tests backbone and terminal alpha_constraints in debug mode, both with and without PPP preprocessing.
    Requires a specific PDB file with two chains, so do not ever change the system here.

    This test is redundant but let's keep it for as long as we do not deprecate CA constraints in PPP completely.
    """
    errors = []
    job = main.run_platform_from_yaml(input_yaml)
    errors = tests.utils.check_file(
        job[0].pele_dir, "pele.conf", terminal_lines_file + default_interval_lines, errors
    )
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
        (7, 0.7, 5.0, interval7_lines, terminal_lines_raw),
        (15, 0.7, 77.0, custom_lines, custom_ter_lines),
    ],
)
def test_alpha_constraints(
    interval, backbone_spring, terminal_spring, expected_ca_lines, expected_ter_lines
):
    """
    Unit test of the alpha carbon constraint builder with various parameters.
    """
    obj = alpha_constraints.AlphaConstraints(
        os.path.join(test_path, "constraints/ca_constraints.pdb"),
        interval=interval,
        backbone_spring=backbone_spring,
        terminal_spring=terminal_spring,
    )
    obj.build_constraints()

    assert obj.terminal_constraints == expected_ter_lines

    for line in expected_ca_lines:
        assert line in obj.backbone_constraints


@pytest.mark.parametrize(
    ("yaml_file", "expected"),
    [
        (
            (None, None, None, 7, "rescoring"),
            set1_lines,
        ),  # should run defaults for rescoring with ca_interval = 7
        (
            (2, 0.7, 7.0, 77, "induced_fit_fast"),
            set2_lines,
        ),  # should run everything with custom values
        ((3, 11.0, None, 11, "rescoring"), set3_lines)
        # should run level 3 with custom ca_interval = 11 and ca_constr spring = 11
    ],
    indirect=["yaml_file"],
)
def test_ca_constraint_logic(yaml_file, expected):
    """
    Tests if the whole EnviroBuilder manages to follow the carbon alpha constraints logic where user arguments >
    constraint levels > package parameters > default values (corresponding to level 1).
    """
    errors = []
    job = main.run_platform_from_yaml(yaml_file)
    errors = tests.utils.check_file(job[0].pele_dir, "pele.conf", expected, errors)

    os.remove(yaml_file)
    assert not errors


@pytest.fixture
def yaml_file(request):
    """
    Unpacks parametrized values for test_ca_constraint_logic and creates a YAML file on the fly for a specific,
    hard-coded system.
    """
    file_name = "input.yaml"
    level, ca_constr, ter_constr, ca_interval, simulation_type = request.param

    args = {
        "system": "../pele_platform/Examples/Msm/PR_1A28_xray_-_minimized.pdb",
        "chain": "Z",
        "resname": "STR",
        simulation_type: True,
        "constraint_level": level,
        "terminal_constr": ter_constr,
        "ca_constr": ca_constr,
        "ca_interval": ca_interval,
        "debug": True,
    }

    with open(file_name, "w+") as f:
        yaml.dump(args, f)

    return file_name


def test_external_constraints(ext_args=EXTERNAL_CONSTR_ARGS):
    errors = []
    job = main.run_platform_from_yaml(ext_args)
    errors = tests.utils.check_file(job[0].pele_dir, "pele.conf", EXT_CONSTR, errors)
    assert not errors


SMILES_CONSTR = [
    '{ "type": "constrainAtomToPosition", "springConstant": 33.5, "equilibriumDistance": 0.0, "constrainThisAtom": "Z:305:_C7_" },',
    '{ "type": "constrainAtomToPosition", "springConstant": 33.5, "equilibriumDistance": 0.0, "constrainThisAtom": "Z:305:_N1_" },',
    '{ "type": "constrainAtomToPosition", "springConstant": 33.5, "equilibriumDistance": 0.0, "constrainThisAtom": "Z:305:_C1_" },',
    '{ "type": "constrainAtomToPosition", "springConstant": 33.5, "equilibriumDistance": 0.0, "constrainThisAtom": "Z:305:_C2_" },',
    '{ "type": "constrainAtomToPosition", "springConstant": 33.5, "equilibriumDistance": 0.0, "constrainThisAtom": "Z:305:_N2_" },',
    '{ "type": "constrainAtomToPosition", "springConstant": 33.5, "equilibriumDistance": 0.0, "constrainThisAtom": "Z:305:_C3_" },',
    '{ "type": "constrainAtomToPosition", "springConstant": 33.5, "equilibriumDistance": 0.0, "constrainThisAtom": "Z:305:_C4_" },',
    '{ "type": "constrainAtomToPosition", "springConstant": 33.5, "equilibriumDistance": 0.0, "constrainThisAtom": "Z:305:_C5_" },',
    '{ "type": "constrainAtomToPosition", "springConstant": 33.5, "equilibriumDistance": 0.0, "constrainThisAtom": "Z:305:_C6_" },',
    '{ "type": "constrainAtomToPosition", "springConstant": 33.5, "equilibriumDistance": 0.0, "constrainThisAtom": "Z:305:_O1_" },',
]


@pytest.mark.parametrize(
    ("yaml_input", "expected"),
    [
        (EXTERNAL_CONSTR_ARGS, EXT_CONSTR),
        (ARGS_SMARTS_CONSTR, SMILES_CONSTR),
    ],
)
def test_constraints(yaml_input, expected):
    """
    Runs platform in debug mode to check if all constraints were correctly set in pele.conf.
    """
    output = main.run_platform_from_yaml(yaml_input)
    errors = []
    folder = output[0].pele_dir if type(output) == list else output.pele_dir
    errors = tests.utils.check_file(folder, "pele.conf", expected, errors)
    assert not errors


@pytest.mark.parametrize(
    "mol_string",
    ["CN1CC[NH+](CC1)CCO", "[#6]-[#7]1-[#6]-[#6]-[#7H+](-[#6]-[#6]-1)-[#6]-[#6]-[#8]"],
)
def test_smiles_constraints_class(mol_string):
    """
    Checks if SmilesConstraints class works correctly, i.e.: converts SMILES to SMARTS, extracts ligand from the PDB,
    matches the right substructure and sets correct constraints.
    """
    obj = smi.SmilesConstraints(
        "../pele_platform/Examples/constraints/4qnr_prep.pdb",
        "CN1CC[NH+](CC1)CCO",
        "LIG",
        "Z",
        33.5,
    )
    smarts_from_smiles = obj.convert_to_smarts(mol_string)


def test_constrain_smarts():
    yaml = os.path.join(test_path, "constraints/input_constrain_smarts.yaml")
    errors = []
    job = main.run_platform_from_yaml(yaml)
    errors = tests.utils.check_file(job[0].pele_dir, "pele.conf", SMILES_CONSTR, errors)
    assert not errors


def test_SmilesConstraints_class():
    obj = smi.SmilesConstraints(
        "../pele_platform/Examples/constraints/4qnr_prep.pdb",
        "CN1CC[NH+](CC1)CCO",
        "LIG",
        "Z",
        33.5,
    )
    smarts_from_smiles = obj.convert_to_smarts("CN1CC[NH+](CC1)CCO")
    smarts_from_smarts = obj.convert_to_smarts(
        "[#6]-[#7]1-[#6]-[#6]-[#7H+](-[#6]-[#6]-1)-[#6]-[#6]-[#8]"
    )
    ligand = obj.extract_ligand(obj.input_pdb, obj.resname)
    matches = obj.get_matches(smarts_from_smiles, ligand)
    constraints = obj.build_constraints(matches, ligand, obj.spring_constant, obj.chain)

    assert (
        smarts_from_smiles == "[#6]-[#7]1-[#6]-[#6]-[#7H+](-[#6]-[#6]-1)-[#6]-[#6]-[#8]"
    )
    assert (
        smarts_from_smiles == "[#6]-[#7]1-[#6]-[#6]-[#7H+](-[#6]-[#6]-1)-[#6]-[#6]-[#8]"
    )
    assert matches == ((9, 0, 1, 2, 3, 4, 5, 6, 7, 8),)
    assert constraints == SMILES_CONSTR
