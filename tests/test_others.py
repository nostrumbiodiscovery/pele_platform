import glob
import os
import pytest
import shutil

import pele_platform.constants.constants as cs
import pele_platform.main as main
import pele_platform.Utilities.Helpers.protein_wizard as pp
import pele_platform.Errors.custom_errors as ce
from pele_platform.Utilities.Helpers import helpers
from pele_platform.Frag import checker as ch

test_path = os.path.join(cs.DIR, "Examples")
EXTERNAL_CONSTR_ARGS = os.path.join(
    test_path, "constraints/input_external_constraints.yaml"
)

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
MULTIPLE_SIMULATIONS = os.path.join(test_path, "checker", "multiple_simulations.yaml")

MAPPED = ['atoms": { "ids":["Z:1:_C13"]}']

EXT_CONSTR = [
    '{ "type": "constrainAtomToPosition", "springConstant": 5, "equilibriumDistance": 0.0, "constrainThisAtom": "A:1:_H__" },',
    '{"type": "constrainAtomsDistance", "springConstant": 50, "equilibriumDistance": 2.34, "constrainThisAtom":  "A:1:_H__", "toThisOtherAtom": "L:1:_C21"}',
]

PPP_CONSTR = [
    '"constrainThisAtom": "B:207:_CA_" }',
    '"constrainThisAtom": "A:1:_CA_" }',
    '"constrainThisAtom": "B:247:_CA_" }',
]


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
        (INPUT_TEMPLATE, ce.TemplateFileNotFound, "File mgadeaz not found"),  # PLOP
        (WRONG_YAML, ce.WrongYamlFile, None),
        (
            UNK_LIGAND,
            ce.LigandNameNotSupported,
            "'UNK' ligand name is not supported, please rename it, e.g. 'LIG'.",
        ),
        (ROTAMER, ce.RotamersFileNotFound, "File mgadeaz not found"),  # PLOP
        (
            OUT_IN,
            ce.OutInError,
            "Flag final_site must be specified for out_in package.",
        ),
        (
            UNDERSCORE,
            ce.WrongAtomStringFormat,
            "Atom string set in A_106:OH does not seem to have the right format. It should follow chain:residue number:atom name pattern, e.g. 'A:105:CA'.",
        ),
        (
            SUBSTRUCTURE_ERROR,
            ce.SubstructureError,
            "More than one substructure found in your ligand. Make sure SMILES constrain pattern is not ambiguous!",
        ),
        (MULTIPLE_SIMULATIONS, ce.MultipleSimulationTypes, None)
    ],
)
def test_errors(yaml_file, expected_error, error_string):
    """
    Checks if all custom errors are raised when expected, e.g. invalid flags in input.yaml or wrong atom names in PDB,
    and the correct error message is printed, if applicable.
    """
    with pytest.raises(expected_error, match=error_string):
        main.run_platform_from_yaml(yaml_file)


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

    with pytest.raises(Exception, match=error):
        main.run_platform_from_yaml(yaml_file)


@pytest.mark.parametrize(
    "complex",
    [
        os.path.join(test_path, "preparation/6qmk_correct.pdb"),
        os.path.join(test_path, "preparation/6qmk_repeated.pdb"),
    ],
)
def test_protein_wizard(complex):
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
    assert len(set(lig_atomnames)) == len(
        lig_atomnames
    )

    assert len(set(repeated_lig_chains)) == 1  # make sure all lig chains are Z
    assert len(set(repeated_lig_atomnames)) == len(repeated_lig_atomnames)


def test_checker_nlimits():

    sdf = os.path.join(test_path, "frag/checker/receptor.sdf")

    with pytest.raises(ce.LigandSizeExceed):
        ch.check_limit_number_atoms(sdf, 100)


def test_checker_subsearch():
    from rdkit import Chem

    ligand = os.path.join(test_path, "frag/checker/4RFM_series.sdf")
    core = os.path.join(test_path, "frag/checker/4RFM_proc.pdb")

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

    with pytest.raises(ce.ExecutableNotInPath):
        main.run_platform_from_yaml(ext_args)

    os.environ["PATH"] = path_variables


def test_skip_preprocess():
    """
    Runs a simulation in debug mode and asserts the skip_preprocess flag works by ensuring the "original" PDB (not
    complex_processed_processed.pdb) is present in the simulation directory.
    """
    yaml = os.path.join(test_path, "box/input_gpcr.yaml")
    (job,) = main.run_platform_from_yaml(yaml)
    input_files = glob.glob(os.path.join(job.inputs_dir, "*_processed.pdb"))
    assert not input_files


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


@pytest.mark.parametrize(
    ("dir_index", "pele_dir"), [(0, "XXX_Pele"), (3, "XXX_Pele_2")]
)
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
