import os
import glob
import pytest
import shutil

import pele_platform.constants.constants as cs
import pele_platform.main as main
from tests import test_adaptive as td

test_path = os.path.join(cs.DIR, "Examples")

FRAG_ARGS = os.path.join(test_path, "frag/input.yaml")
FRAG_SIM_ARGS = os.path.join(test_path, "frag/input_sim.yaml")
FRAG_CORE_ARGS = os.path.join(test_path, "frag/input_core.yaml")
FLAGS_ARGS = os.path.join(test_path, "frag/input_flags.yaml")
FRAG_JOINER_ARGS = os.path.join(test_path, "frag/sdf_joiner/*.yml")
FRAG_SDF_LIBRARIES = os.path.join(test_path, "frag/input_lib_sdf.yaml")
FRAG_PDB_LIBRARIES = os.path.join(test_path, "frag/input_lib_pdb.yaml")
FRAG_ANALYSIS_TO_POINT = os.path.join(test_path, "frag/input_point_analysis.yaml")
FRAG_SYMMETRY = os.path.join(test_path, "frag/input_symmetry.yaml")

EXPECTED_INPUT = os.path.join(
    test_path, "frag/asymmetric_hydrogens_detector/expected_input.conf"
)

PDB_lines = [
    "../pele_platform/Examples/frag/pdb_test_lib/mol1.pdb C3-H2 C1-H1",
    "../pele_platform/Examples/frag/pdb_test_lib/mol1.pdb C3-H2 C2-H4",
    "../pele_platform/Examples/frag/pdb_test_lib/mol1.pdb C3-H2 N1-H6",
    "../pele_platform/Examples/frag/pdb_test_lib/mol2.pdb C3-H2 C1-H1",
    "../pele_platform/Examples/frag/pdb_test_lib/mol2.pdb C3-H2 N1-H4",
]

SDF_lines = [
    "test_library-1.pdb C3-H2 C1-H1",
    "test_library-1.pdb C3-H2 C2-H4",
    "test_library-1.pdb C3-H2 N1-H6",
    "test_library-2.pdb C3-H2 C1-H1",
    "test_library-2.pdb C3-H2 N1-H4",
]

point_analysis_lines = [
    "../pele_platform/Examples/frag/analysis_data/1w7h_preparation_structure_2w_processed_mol1C3-H2C1-H1,../pele_platform/Examples/frag/analysis_data/1w7h_preparation_structure_2w_processed_mol1C3-H2C1-H1/top_result/epochsampling_result_trajectory_1.1_BindingEnergy-23.4636.pdb,2.73029273852075,../pele_platform/Examples/frag/analysis_data/1w7h_preparation_structure_2w_processed_mol1C3-H2C1-H1/top_result/epochsampling_result_trajectory_2.1_BindingEnergy-25.1634.pdb,-25.1634,../pele_platform/Examples/frag/analysis_data/1w7h_preparation_structure_2w_processed_mol1C3-H2C1-H1/top_result/epochsampling_result_trajectory_1.1_BindingEnergy-23.4636.pdb,0.7087330424726306,2.730292738520753,-23.4636"
]


def test_frag_sim(
    capsys,
    ext_args=FRAG_SIM_ARGS,
    output="1w7h_preparation_structure_2w_processed_aminoC1N1",
):
    """
    Runs FragPELE test simulation. Checks if the output folder exists and the ligand was not skipped.

    Parameters
    ----------
    ext_args : str
        Path to PELE input file.
    output : str
        Output folder name.
    """

    if os.path.exists(output):
        shutil.rmtree(output)

    job = main.run_platform(ext_args)
    captured = capsys.readouterr()
    top_results = glob.glob(os.path.join(output, "top_result", "*pdb"))

    assert "Skipped - FragPELE will not run." not in captured.out
    assert os.path.exists(output)
    assert len(top_results) == 3


def test_frag_core(capsys, ext_args=FRAG_CORE_ARGS):
    """
    Tests FragPELE growing method using an SDF with full ligands. Checks if the output folder exists and the ligand
    was not skipped.

    Parameters
    ----------
    ext_args : str
        Path to PELE input file.
    """

    output = "1w7h_preparation_structure_2w_processed_*ligand_01*"

    output_paths = glob.glob(output)
    for path in output_paths:
        if os.path.exists(output):
            shutil.rmtree(output)

    job = main.run_platform(ext_args)
    captured = capsys.readouterr()

    new_output_path = glob.glob(output)[0]
    top_results = glob.glob(os.path.join(new_output_path, "top_result", "*pdb"))

    assert "Skipped - FragPELE will not run." not in captured.out
    assert os.path.exists(new_output_path)
    assert len(top_results) == 3


def test_flags(ext_args=FLAGS_ARGS, output="water_processed_aminoCA1N1"):
    """
    Checks input file flags.

    Parameters
    ----------
    ext_args : str
        Path to PELE input file.
    output : str
        Output folder name.
    """
    FRAG_FLAGS = [
        '"seed" : 3000',
    ]
    errors = []
    if os.path.exists(output):
        shutil.rmtree(output, ignore_errors=True)
    job = main.run_platform(ext_args)
    folder = output
    errors = td.check_file(
        folder,
        "control_folder/0_pele_template.conf",
        td.PELE_VALUES + FRAG_FLAGS,
        errors,
    )
    errors = td.check_file(
        folder, "DataLocal/LigandRotamerLibs/SB4.rot.assign", "60", errors
    )
    assert not errors


def test_sdf_joiner(ext_args=FRAG_JOINER_ARGS):
    """
    Tests the SDF joiner.

    Parameters
    ----------
    ext_args : str
        Path to PELE input file.
    """
    files = glob.glob(ext_args)
    for file in files:
        try:
            job = main.run_platform(file)
        except Exception:
            assert False


@pytest.mark.parametrize(
    ("yaml_file", "expected_lines"),
    [(FRAG_SDF_LIBRARIES, SDF_lines), (FRAG_PDB_LIBRARIES, PDB_lines)],
)
def test_libraries(yaml_file, expected_lines):
    """
    Tests the growing of fragments from a custom-made SDF and PDB libraries.

    Parameters
    ----------
    yaml_file : str
        Path to PELE input file.
    expected_lines : list[str]
        List of lines expected in input.conf.
    """
    if os.path.exists("input.conf"):
        os.remove("input.conf")

    job = main.run_platform(yaml_file)
    errors = []
    errors = td.check_file(os.getcwd(), "input.conf", expected_lines, errors)
    assert not errors


def test_analysis_to_point(ext_args=FRAG_ANALYSIS_TO_POINT):
    """
    Tests the automated analysis to retrieve most promising fragments
    from a custom-made library based on their proximity to a certain point.

    Parameters
    ----------
    ext_args : str
        Path to PELE input file.
    """
    job = main.run_platform(ext_args)
    errors = []
    errors = td.check_file(
        os.getcwd(), "point_analysis.csv", point_analysis_lines, errors, ",", 4
    )
    assert not errors


def test_symmetry(ext_args=FRAG_SYMMETRY):
    """
    Tests the asymmetric hydrogen detector.

    Parameters
    ----------
    ext_args : str
        Path to PELE input file.
    """
    if os.path.exists("input.conf"):
        os.remove("input.conf")
    job = main.run_platform(ext_args)
    errors = []
    errors = td.check_file(os.getcwd(), "input.conf", PDB_lines, errors)
    assert not errors
