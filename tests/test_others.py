import os
import pytest
from subprocess import Popen, PIPE

import pele_platform.constants.constants as cs
import pele_platform.main as main
import pele_platform.Utilities.Helpers.protein_wizard as pp
import pele_platform.Errors.custom_errors as ce
from pele_platform.Utilities.Helpers import smiles_constraints as smi
from . import test_adaptive as tk

test_path = os.path.join(cs.DIR, "Examples")

EXTERNAL_CONSTR_ARGS = os.path.join(test_path, "constraints/input_external_constraints.yaml")
PPP_CONSTR_ARGS = os.path.join(test_path, "constraints/input_ppp.yaml")
ARGS_SMARTS_CONSTR = os.path.join(test_path, "constraints/input_constrain_smarts.yaml")
LIG_PREP_ARGS = os.path.join(test_path, "preparation/input_space.yaml")
ENV_ARGS = os.path.join(test_path, "checker/input_env.yaml")
ATOM_GPCR_ERROR_ARGS = os.path.join(test_path, "gpcr/input_atom_error.yaml")
LIGAND_CHECKER = os.path.join(test_path, "checker/input.yaml")
SDF_LIMIT_ATOMS = os.path.join(test_path, "frag/checker/receptor.sdf")
ATOM_STRING = os.path.join(test_path, "checker/input_atom_string.yaml")
FLAG_SIMILARITY = os.path.join(test_path, "checker/input.yaml")
INPUT_TEMPLATE = os.path.join(test_path, "checker/input_template.yaml")
WRONG_YAML = os.path.join(test_path, "gpcr/complex.pdb")
UNK_LIGAND = os.path.join(test_path, "checker/input_unk.yaml")
ROTAMER = os.path.join(test_path, "checker/input_rotamer.yaml")
OUT_IN = os.path.join(test_path, "out_in/input_flag_error.yaml")
UNDERSCORE = os.path.join(test_path, "checker/input_underscore.yaml")
SUBSTRUCTURE_ERROR = os.path.join(test_path, "constraints/input_smiles_error.yaml")

EXT_CONSTR = [
    '{ "type": "constrainAtomToPosition", "springConstant": 5, "equilibriumDistance": 0.0, "constrainThisAtom": "A:1:_H__" },',
    '{"type": "constrainAtomsDistance", "springConstant": 50, "equilibriumDistance": 2.34, "constrainThisAtom":  "A:1:_H__", "toThisOtherAtom": "L:1:_C21"}',
]

PPP_CONSTR = [
    '"constrainThisAtom": "B:207:_CA_" }',
    '"constrainThisAtom": "A:1:_CA_" }',
    '"constrainThisAtom": "B:247:_CA_" }',
]

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
    ("yaml_file", "expected"),
    [
        (EXTERNAL_CONSTR_ARGS, EXT_CONSTR),
        (PPP_CONSTR_ARGS, PPP_CONSTR),
        (ARGS_SMARTS_CONSTR, SMILES_CONSTR),
    ],
)
def test_constraints(yaml_file, expected):
    """
    Runs platform in debug mode to check if all constraints were correctly set in pele.conf.
    """
    output = main.run_platform(yaml_file)
    errors = []
    folder = output[0].pele_dir if type(output) == list else output.pele_dir
    errors = tk.check_file(folder, "pele.conf", expected, errors)
    assert not errors


@pytest.mark.parametrize(
    ("yaml_file", "expected_error", "error_string"),
    [
        (LIGAND_CHECKER, KeyError, None),
        (LIG_PREP_ARGS, ce.LigandPreparationError, None),
        (ENV_ARGS, ce.EnvVariableNotFound, None),
        (ATOM_STRING, ce.WrongAtomStringFormat, None),
        (FLAG_SIMILARITY, KeyError, "Incorrect flag posis. Did you mean poses?"),
        (
                ATOM_GPCR_ERROR_ARGS,
                ce.WrongAtomSpecified,
                "Atom A:114:CM could not be found in structure",
        ),
        (INPUT_TEMPLATE, ce.TemplateFileNotFound, "File mgadeaz not found"),
        (WRONG_YAML, ce.WrongYamlFile, None),
        (
                UNK_LIGAND,
                ce.LigandNameNotSupported,
                "'UNK' ligand name is not supported, please rename it, e.g. 'LIG'.",
        ),
        (ROTAMER, ce.RotamersFileNotFound, "File mgadeaz not found"),
        (
                OUT_IN,
                ce.OutInError,
                "Flag final_site must be specified for out_in package.",
        ),
        (
                UNDERSCORE,
                ce.WrongAtomStringFormat,
                "The specified atom is wrong 'A_106:OH'. Should be 'chain:resnumber:atomname",
        ),
        (
                SUBSTRUCTURE_ERROR,
                ce.SubstructureError,
                "More than one substructure found in your ligand. Make sure SMILES constrain pattern is not ambiguous!",
        ),
    ],
)
def test_errors(yaml_file, expected_error, error_string):
    """
    Checks if all custom errors are raised when expected, e.g. invalid flags in input.yaml or wrong atom names in PDB,
    and the correct error message is printed, if applicable.
    """
    try:
        job = main.run_platform(yaml_file)
    except expected_error as e:
        assert expected_error
        if error_string:
            assert str(e) == error_string  # or str(e).strip("'") == error_string  # stupid, but works
        return
    assert False


@pytest.mark.parametrize("complex", [os.path.join(test_path, "preparation/6qmk_correct.pdb"),
                                     os.path.join(test_path, "preparation/6qmk_repeated.pdb")])
def test_proteinwizard(complex):
    """
    Checks if protein wizard function correctly preprocesses the protein, sets a unique ligand chain, etc.
    """
    yaml = os.path.join(test_path, "preparation/input.yaml")
    output = pp.prep_complex(complex, yaml, prep_output=complex, debug=True)

    with open(output, "r") as f:
        lines = f.readlines()
        for line in lines:
            lig_chains = [line[21:22].strip() for line in lines if line[17:20] == "J8H"]
            lig_atomnames = [line[12:16].strip() for line in lines if line[17:20] == "J8H"]

    assert all(elem == "Z" for elem in lig_chains)
    assert len(set(lig_atomnames)) == len(lig_atomnames)


def test_mpirun_in_path(ext_args=EXTERNAL_CONSTR_ARGS):
    path_variables = os.environ["PATH"]
    os.environ["PATH"] = ""
    try:
        job = main.run_platform(ext_args)
    except ce.ExecutableNotInPath:
        assert True
        os.environ["PATH"] = path_variables
        return
    os.environ["PATH"] = path_variables
    assert False


def test_python_version_error():
    p = Popen("/usr/bin/python -m pele_platform.main -h".split(), stdout=PIPE, stderr=PIPE)
    output, error = p.communicate()
    if "OldPythonVersion" in error.decode():
        assert True
        return
    assert False


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
    ligand = obj.extract_ligand(obj.input_pdb, obj.resname)
    matches = obj.get_matches(smarts_from_smiles, ligand)
    constraints = obj.build_constraints(matches, ligand, obj.spring_constant, obj.chain)

    assert smarts_from_smiles == "[#6]-[#7]1-[#6]-[#6]-[#7H+](-[#6]-[#6]-1)-[#6]-[#6]-[#8]"
    assert matches == ((9, 0, 1, 2, 3, 4, 5, 6, 7, 8),)
    assert constraints == SMILES_CONSTR


def test_skip_preprocess():
    """
    Runs a simulation in debug mode and asserts the skip_preprocess flag works by ensuring the "original" PDB (not
    complex_processed_processed.pdb) is present in the simulation directory.
    """
    yaml = os.path.join(test_path, "box/input_gpcr.yaml")
    (job,) = main.run_platform(yaml)
    input_file = os.path.join(job.pele_dir, "complex_processed.pdb")
    assert os.path.exists(input_file)
