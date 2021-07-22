import os
import pytest
from subprocess import Popen, PIPE

import shutil

import pele_platform.constants.constants as cs
import pele_platform.main as main
import pele_platform.Utilities.Helpers.protein_wizard as pp
import pele_platform.Errors.custom_errors as ce
from pele_platform.Utilities.Helpers.constraints import smiles_constraints as smi
from pele_platform.Utilities.Helpers import helpers

from pele_platform.Utilities.Helpers import smiles_constraints as smi
from . import test_adaptive as tk

test_path = os.path.join(cs.DIR, "Examples")
EXTERNAL_CONSTR_ARGS = os.path.join(
    test_path, "constraints/input_external_constraints.yaml"
)

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
MAP_ARGS = os.path.join(test_path, "checker/input_map_atom_str.yaml")

MAPPED = ['atoms": { "ids":["Z:1:_C13"]}']

EXT_CONSTR = [
    '{ "type": "constrainAtomToPosition", "springConstant": 5, "equilibriumDistance": 0.0, "constrainThisAtom": "A:1:_H__" },',
    '{"type": "constrainAtomsDistance", "springConstant": 50, "equilibriumDistance": 2.34, "constrainThisAtom":  "A:1:_H__", "toThisOtherAtom": "L:1:_C21"}',
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


def test_external_constraints(ext_args=EXTERNAL_CONSTR_ARGS):
    errors = []
    job = main.run_platform_from_yaml(ext_args)
    errors = tk.check_file(job.pele_dir, "pele.conf", EXT_CONSTR, errors)
    assert not errors

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
            "The specified atom is wrong 'A_106:OH'. Should be 'chain:resnumber:atomname'",
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
        job = main.run_platform_from_yaml(yaml_file)
    except expected_error as e:
        assert expected_error
        if error_string:
            assert (str(e) == error_string) or str(e).strip(
                "'"
            ) == error_string  # stupid, but works
        return
    assert False


@pytest.mark.parametrize(
    ("yaml_file", "error"),
    [
        (
            os.path.join(test_path, "checker/input_space.yaml"),
            "Please ensure every key in input.yaml is followed by a colon and a space. There seem to be some issues on line 2, character 1.",
        ),
        (
            os.path.join(test_path, "checker/input_tab.yaml"),
            "Please remove any trailing tabs from input.yaml, there seem to be one on line 1, character 60.",
        ),
    ],
)
def test_yaml_errors(yaml_file, error):
    try:
        job = main.run_platform_from_yaml(yaml_file)
    except Exception as e:
        assert str(e) == error
        return
    assert False


def test_proteinwizard():
    complex_correct = os.path.join(test_path, "preparation/6qmk_correct.pdb")
    complex_repeated = os.path.join(test_path, "preparation/6qmk_repeated.pdb")
@pytest.mark.parametrize(
    "complex",
    [
        os.path.join(test_path, "preparation/6qmk_correct.pdb"),
        os.path.join(test_path, "preparation/6qmk_repeated.pdb"),
    ],
)
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
            lig_atomnames = [
                line[12:16].strip() for line in lines if line[17:20] == "J8H"
            ]
            repeated_lig_chains = [
                line[21:22].strip() for line in lines if line[17:20] == "J8H"
            ]
            repeated_lig_atomnames = [
                line[12:16].strip() for line in lines if line[17:20] == "J8H"
            ]

    assert all(elem == "Z" for elem in lig_chains)
    assert len(set(lig_atomnames)) == len(lig_atomnames)  # make sure all lig chains are Z
    assert len(set(repeated_lig_atomnames)) == len(repeated_lig_atomnames)


SDF_LIMIT_ATOMS = os.path.join(test_path, "frag/checker/receptor.sdf")


def test_checker_nlimits(sdf=SDF_LIMIT_ATOMS):
    try:
        ch.check_limit_number_atoms(sdf, 100)
    except ce.LigandSizeExceed:
        assert True
    else:
        assert False


PDB_CHECKER_SUBSEARCH = os.path.join(test_path, "frag/checker/4RFM_series.sdf")
CORE_CHECKER_SUBSEARCH = os.path.join(test_path, "frag/checker/4RFM_proc.pdb")


def test_checker_subsearch(ligand=PDB_CHECKER_SUBSEARCH, core=CORE_CHECKER_SUBSEARCH):
    from rdkit import Chem

    mol = next(Chem.SDMolSupplier(ligand))
    core = Chem.rdmolfiles.MolFromPDBFile(core)
    atoms_in_common = mol.GetSubstructMatches(core)[0]
    atoms_in_common_after = ch.check_substructure_match(core, mol, atoms_in_common)
    # Exchange nitrogen due to wrong previous result
    assert atoms_in_common != atoms_in_common_after
    assert atoms_in_common_after[atoms_in_common.index(13)] == 12
    assert atoms_in_common_after[atoms_in_common.index(12)] == 13


def test_mpirun_in_path(ext_args=EXTERNAL_CONSTR_ARGS):
    path_variables = os.environ["PATH"]
    os.environ["PATH"] = ""
    try:
        job = main.run_platform_from_yaml(ext_args)
    except ce.ExecutableNotInPath:
        assert True
        os.environ["PATH"] = path_variables
        return
    os.environ["PATH"] = path_variables
    assert False


def test_lig_preparation_error(args=LIG_PREP_ARGS):
    try:
        job = main.run_platform_from_yaml(args)
    except ce.LigandPreparationError:
        assert True
        return
    assert False


def test_env_variable(ext_args=ENV_ARGS):
    try:
        job = main.run_platform_from_yaml(ext_args)
    except ce.EnvVariableNotFound as e:
        assert True
        return
    assert False



def test_flag_similarity():
    yaml = os.path.join(test_path, "checker/input.yaml")
    try:
        job = main.run_platform_from_yaml(yaml)
    except KeyError as e:
        assert str(e).strip("'") == "Incorrect flag posis. Did you mean poses?"
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


def test_atom_error(ext_args=ATOM_GPCR_ERROR_ARGS):
    try:
        job = main.run_platform_from_yaml(ext_args)
    except ce.WrongAtomSpecified as e:
        assert str(e).strip("'") == "Atom A:114:CM could not be found in structure"
        return
    assert False


yaml = os.path.join(test_path, "checker/input_template.yaml")


def test_template_error(yaml=yaml):
    try:
        job = main.run_platform_from_yaml(yaml)
    except ce.TemplateFileNotFound as e:
        assert (
            str(e).strip("'")
            == "Could not locate mgadeaz file. Please double-check the path."
        )
        return
    assert False


def test_input_yaml_error():
    yaml = os.path.join(test_path, "gpcr/complex.pdb")
    try:
        job = main.run_platform_from_yaml(yaml)
    except ce.WrongYamlFile:
        assert True
        return
    assert False


yaml = os.path.join(test_path, "checker/input_rotamer.yaml")


def test_rotamer_error(yaml=yaml):
    try:
        job = main.run_platform_from_yaml(yaml)
    except ce.RotamersFileNotFound as e:
        assert (
            str(e).strip("'")
            == "Could not locate mgadeaz file. Please double-check the path."
        )
        return
    assert False


yaml = os.path.join(test_path, "out_in/input_flag_error.yaml")


def test_out_in_flag(yaml=yaml):
    try:
        job = main.run_platform_from_yaml(yaml)
    except ce.OutInError as e:
        assert (
            str(e).strip("'") == "flag final_site must be specified for out_in package"
        )
        return
    assert False


yaml = os.path.join(test_path, "checker/input_atom_string.yaml")


def test_atom_string_error(yaml=yaml):
    try:
        job = main.run_platform_from_yaml(yaml)
    except ce.WrongAtomStringFormat as e:
        assert (
            str(e).strip("'")
            == "The specified atom is wrong '157:A:N'. Should be 'chain:resnumber:atomname"
        )
        return
    assert False


yaml = os.path.join(test_path, "checker/input_underscore.yaml")


def test_atom_string_underscore(yaml=yaml):
    try:
        job = main.run_platform_from_yaml(yaml)
    except ce.WrongAtomStringFormat as e:
        assert (
            str(e).strip("'")
            == "The specified atom is wrong 'A_106:OH'. Should be 'chain:resnumber:atomname"
        )


yaml = os.path.join(test_path, "checker/input_unk.yaml")


def test_unk_error():
    try:
        job = main.run_platform_from_yaml(yaml)
    except ce.LigandNameNotSupported as e:
        assert (
            str(e)
            == "'UNK' ligand name is not supported, please rename it, e.g. 'LIG'."
        )
        return
    assert False


def test_constrain_smarts():
    yaml = os.path.join(test_path, "constraints/input_constrain_smarts.yaml")
    errors = []
    job = main.run_platform_from_yaml(yaml)
    errors = tk.check_file(job.pele_dir, "pele.conf", SMILES_CONSTR, errors)
    assert not errors


def test_substructure_error():
    yaml = os.path.join(test_path, "constraints/input_smiles_error.yaml")
    try:
        job = main.run_platform_from_yaml(yaml)
    except ce.SubstructureError as e:
        assert (
            str(e).strip("'")
            == "More than one substructure found in your ligand. Make sure SMILES constrain pattern is not ambiguous!"
        )


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
        smarts_from_smiles
        == "[#6]-[#7]1-[#6]-[#6]-[#7H+](-[#6]-[#6]-1)-[#6]-[#6]-[#8]"
    )
    assert (
        smarts_from_smiles == "[#6]-[#7]1-[#6]-[#6]-[#7H+](-[#6]-[#6]-1)-[#6]-[#6]-[#8]"
    )
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


def test_check_multiple_simulations():
    """
    Ensures the platform raises an error, if the user sets more than one simulation type in YAML.
    """
    yaml_file = os.path.join(test_path, "checker", "multiple_simulations.yaml")

    with pytest.raises(ce.MultipleSimulationTypes):
        main.run_platform_from_yaml(yaml_file)


@pytest.mark.parametrize(
    ("pdb_file", "residues", "expected"),
    [
        ("1zop.pdb", ["MG"], {"MG": [" MG "]}),
        ("4w4o_prep.pdb", ["ZN"], {"ZN": ["ZN  "]}),
    ],
)
def test_retrieve_atom_names(pdb_file, residues, expected):
    """
    Tests extraction of PDB atom names from a PDB file for a specific set of residues.

    Parameters
    -----------
    pdb_file : str
        File name in Examples/constraints.
    residues : List[str]
        List of residue names to extract.
    expected : dict
        Expected output dictionary, where keys are residue names and values - a list of PDB atom names.
    """
    pdb_file = os.path.join(test_path, "constraints", pdb_file)
    output = helpers.retrieve_atom_names(pdb_file, residues)

    assert output == expected


@pytest.mark.parametrize(("dir_index", "pele_dir"), [(0, "XXX_Pele"), (3, "XXX_Pele_2")])
def test_latest_pele_dir(dir_index, pele_dir):
    """
    Tests if the platform correctly retrieves the latest pele_dir.

    Parameters
    -----------
    dir_index : int
        Index of the expected latest directory.
    pele_dir : str
        PELE directory that would normally be created by the platform.
    """
    # Make folders to mock existing simulation output (mind that range indexes from 0)
    for i in range(dir_index + 1):
        os.mkdir(f"XXX_Pele_{i}")

    output = helpers.get_latest_peledir(pele_dir)
    expected_output = f"XXX_Pele_{dir_index}" if dir_index > 0 else "XXX_Pele"
    assert output == expected_output

    # Clean up
    for i in range(dir_index + 1):
        shutil.rmtree(f"XXX_Pele_{i}")


def test_flag_compatibility_checker():
    """
    Checks if the platform raises an error when incompatible flags are introduced.
    """

    yaml = os.path.join(test_path, "checker", "incompatible_flags.yaml")

    with pytest.raises(ce.IncompatibleYamlFlags):
        main.run_platform_from_yaml(yaml)
