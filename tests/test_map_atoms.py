import os
import pytest

from pele_platform.constants import constants
from pele_platform import main
from . import test_adaptive as tk
from pele_platform.Utilities.Helpers import map_atoms

test_path = os.path.join(constants.DIR, "Examples")

MAP_STR_YAML = os.path.join(test_path, "checker/input_map_atom_str.yaml")
MAP_NUM_YAML = os.path.join(test_path, "checker/input_map_atom_num.yaml")
PDB = os.path.join(test_path, "checker/mapping_complex.pdb")
MAPPED_STR = ['atoms": { "ids":["Z:1:_C13"]}']


@pytest.mark.parametrize(
    ("yaml_file", "expected"), [(MAP_STR_YAML, MAPPED_STR), (MAP_NUM_YAML, MAPPED_STR)]
)
def test_atom_string_mapping(yaml_file, expected):
    """
    Integration tests to ensure both atom strings and atom numbers are correctly mapped and injected into
    pele.conf JSON."
    """
    errors = []
    job = main.run_platform(yaml_file)
    errors = tk.check_file(job.pele_dir, "pele.conf", expected, errors)
    assert not errors


@pytest.mark.parametrize(
    ("line", "expected_atomname", "expected_resnum", "expected_residue_name", "expected_chain"),
    [
        (
            "ATOM     13  HB3 ASP A  11      34.775  63.714  47.018  1.00  0.00           H",
            "HB3",
            "11",
            "ASP",
            "A",
        ),
        (
            "HETATM 5207  N2  API Z 900      -0.107  -0.083   0.156  1.00  0.00           N1+",
            "N2",
            "900",
            "API",
            "Z",
        ),
    ],
)
def test_get_atom_from_line(line, expected_atomname, expected_resnum, expected_residue_name, expected_chain):
    """
    Test to check if extraction of PDB atom name, residue number and chain ID from a PDB line works.
    """
    atom_name, residue_number, residue_name, chain_id = map_atoms.get_atom_from_line(line)
    assert atom_name == expected_atomname
    assert residue_number == expected_resnum
    assert residue_name == expected_residue_name
    assert chain_id == expected_chain


@pytest.mark.parametrize(
    ("line", "expected_output"),
    [
        (
            "ATOM     13  HB3 ASP A  11      34.775  63.714  47.018  1.00  0.00           H",
            ["34.775", "63.714", "47.018"],
        ),
        (
            "HETATM 5207  N2  API Z 900      -0.107  -0.083   0.156  1.00  0.00           N1+",
            ["-0.107", "-0.083", "0.156"],
        ),
    ],
)
def test_get_coords_from_line(line, expected_output):
    """
    Test to check if extraction of atom coordinates from a PDB line works.
    """
    output = map_atoms.get_coords_from_line(line)
    assert output == expected_output


@pytest.mark.parametrize(
    ("atom_number", "expected_output"),
    [
        (428, ["A:40:CA"]),
        (7715, ["Z:201:C1"]),
        ("A:40:CA", ["A:40:CA"]),
        ([7715, "Z:201:C1"], ["Z:201:C1", "Z:201:C1"]),
    ],
)
def test_atom_number_to_atom_string(atom_number, expected_output):

    output = map_atoms.atom_number_to_atom_string(PDB, atom_number)
    assert output == expected_output
