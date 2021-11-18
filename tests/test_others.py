from Bio.PDB import PDBParser
import numpy as np
import os
import shutil

import pele_platform.constants.constants as cs
import pele_platform.constants.pele_params as pp
import pele_platform.main as main
from . import test_adaptive as tk
import pytest
import pele_platform.Utilities.Helpers.protein_wizard as pp
import pele_platform.Frag.checker as ch
import pele_platform.Errors.custom_errors as ce
from pele_platform.Utilities.Helpers.constraints import smiles_constraints as smi
from pele_platform.Utilities.Helpers import helpers
from pele_platform.Utilities.Helpers import randomize

test_path = os.path.join(cs.DIR, "Examples")
EXTERNAL_CONSTR_ARGS = os.path.join(
    test_path, "constraints/input_external_constraints.yaml"
)
LIG_PREP_ARGS = os.path.join(test_path, "preparation/input_space.yaml")
ENV_ARGS = os.path.join(test_path, "checker/input_env.yaml")
ATOM_GPCR_ERROR_ARGS = os.path.join(test_path, "gpcr/input_atom_error.yaml")
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


def test_checker():
    yaml = os.path.join(test_path, "checker/input.yaml")
    try:
        job = main.run_platform_from_yaml(yaml)
    except KeyError:
        assert KeyError
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
    yaml = os.path.join(test_path, "preparation/input.yaml")

    correct_output = pp.prep_complex(
        complex_correct, yaml, prep_output=complex_correct, debug=True
    )
    repeated_output = pp.prep_complex(
        complex_repeated, yaml, prep_output=complex_repeated, debug=True
    )

    with open(correct_output, "r") as f:
        lines = f.readlines()
        for line in lines:
            correct_lig_chains = [
                line[21:22].strip() for line in lines if line[17:20] == "J8H"
            ]
            correct_lig_atomnames = [
                line[12:16].strip() for line in lines if line[17:20] == "J8H"
            ]

    with open(repeated_output, "r") as f:
        lines = f.readlines()
        for line in lines:
            repeated_lig_chains = [
                line[21:22].strip() for line in lines if line[17:20] == "J8H"
            ]
            repeated_lig_atomnames = [
                line[12:16].strip() for line in lines if line[17:20] == "J8H"
            ]

    assert all(elem == "Z" for elem in repeated_lig_chains) and all(
        elem == "Z" for elem in correct_lig_chains
    )  # make sure all lig chains are Z
    assert len(set(correct_lig_atomnames)) == len(
        correct_lig_atomnames
    )  # check for repetition
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
            str(e).strip("'")
            == "flag final_site must be specified for the out_in package."
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
            == "The specified atom is wrong 'A_106:OH'. Should be 'chain:residue number:atom name'."
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
        == smarts_from_smarts
        == "[#6]-[#7]1-[#6]-[#6]-[#7H+](-[#6]-[#6]-1)-[#6]-[#6]-[#8]"
    )
    assert matches == ((9, 0, 1, 2, 3, 4, 5, 6, 7, 8),)
    assert constraints == SMILES_CONSTR


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


def test_flag_compatibility_checker():
    """
    Checks if the platform raises an error when incompatible flags are introduced.
    """

    yaml = os.path.join(test_path, "checker", "incompatible_flags.yaml")

    with pytest.raises(ce.IncompatibleYamlFlags):
        main.run_platform_from_yaml(yaml)


def test_out_in_metrics():
    """
    Production run to test user-defined metrics in OutIn as well as setting the right atom distance for
    the bias column.
    """
    expected_metrics = [
        '"atoms": { "ids":["A:867:_CB_"]}',  # final site
        '"links": { "ids":["Z:1"]}',  # ligand
        '"atoms": { "ids":["A:583:_NH2"]}',  # atom_dist set with atom string
        '"atoms": { "ids":["A:580:_HA_"]}',  # atom_dist set with atom number
    ]
    yaml_file = os.path.join(test_path, "out_in", "input_metrics.yaml")
    job = main.run_platform_from_yaml(yaml_file)
    errors = tk.check_file(job.pele_dir, "pele.conf", expected_metrics, [])
    assert not errors


@pytest.mark.parametrize(
    ("site_finder", "ppi", "center_of_interface", "box_center", "box_radius", "expected_sphere_cent",),
    (
        # classic site_finder: sphere center at protein COM, radius depends on protein size
        [True, False, False, None, None, [58.46355848210691, 7.33174024600783, -24.401180347171945]],
        # site finder focused on a box: sphere cent = box center, radius = box_radius
        [True, False, False, "A:686:CG2", 20, [47.18299865722656, -21.364999771118164, -30.69099998474121]],
        # # classic PPI: sphere cent = COI, default radius = 10
        [False, True, "A:979:OD1", None, None, [53.71099853515625, -28.28499984741211, -8.913999557495117]],
        # # other cases (GPCR, OutIn) without box: sphere cent = COI/initial_site, default radius
        [False, False, "A:979:OD1", None, None, [53.71099853515625, -28.28499984741211, -8.913999557495117]],
        # # other cases (GPCR, OutIn) with box: sphere cent = COI/initial_site, radius = box radius
        [False, False, "A:979:OD1", "A:686:CG2", 5, [53.71099853515625, -28.28499984741211, -8.913999557495117]],
    ),
)
def test_randomize(
    mock_parameters, site_finder, ppi, center_of_interface, box_center, box_radius, expected_sphere_cent,
):
    """
    Checks if randomization works correctly in all possible scenarios.

    Parameters
    ----------
    mock_parameters : Mock
        Mock object to replace Parameters generated during the simulation.
    site_finder : bool
        Whether it's running SiteFinder.
    ppi : bool
        Whether it's running PPI.
    center_of_interface : str
        Center of interface string, e.g. "A:123:CA"
    box_center : str
        Atom string for box center, e.g. "A:123:CA"
    box_radius : float
        Radius of the simulation box.
    """
    mock_parameters.site_finder = site_finder
    mock_parameters.ppi = ppi
    mock_parameters.center_of_interface = center_of_interface

    # Biopython complains if the inputs directory doesn't exist
    if not os.path.exists(mock_parameters.inputs_dir):
        os.mkdir(mock_parameters.inputs_dir)

    output_files, sphere_radius, sphere_center = randomize.randomize_starting_position(
        mock_parameters, box_center=box_center, box_radius=box_radius
    )

    assert sphere_center == expected_sphere_cent

    if box_radius:
        assert sphere_radius >= box_radius

    distances = []
    for file in output_files:
        parser = PDBParser()
        structure = parser.get_structure('protein', file)
        com = randomize.calculate_com(structure)
        distance = np.linalg.norm(com - sphere_center)
        distances.append(distance)

    if box_radius:
        assert all([distance <= box_radius + 1 for distance in distances])
    else:
        assert all([distance <= sphere_radius for distance in distances])
    print("Max distance from sphere center:", max(distances))

    # Clean up after the test
    if os.path.exists(mock_parameters.inputs_dir):
        shutil.rmtree(mock_parameters.inputs_dir)


@pytest.fixture
def mock_parameters():
    folder = os.path.join(test_path, "randomize")
    inputs_dir = "randomize"

    from pele_platform.Utilities.Parameters import Parameters
    import logging

    parameters = Parameters
    parameters.system = os.path.join("site_finder", "complex.pdb")
    parameters.residue = "API"
    parameters.chain = "Z"
    parameters.receptor = os.path.join(folder, "receptor.pdb")
    parameters.ligand_ref = os.path.join(folder, "ligand.pdb")
    parameters.poses = 10
    parameters.inputs_dir = inputs_dir
    parameters.test = True
    parameters.logger = logging.getLogger(__name__)

    return parameters
