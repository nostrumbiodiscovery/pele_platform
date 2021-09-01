import os
import glob
import pytest
import shutil

import pele_platform.constants.constants as cs
import pele_platform.main as main
import tests.utils
from tests import test_adaptive as td

test_path = os.path.join(cs.DIR, "Examples")

FRAG_ARGS = os.path.join(test_path, "frag/input.yaml")
FRAG_CORE_ARGS = os.path.join(test_path, "frag/input_core.yaml")
FLAGS_ARGS = os.path.join(test_path, "frag/input_flags.yaml")
FRAG_JOINER_ARGS = os.path.join(test_path, "frag/sdf_joiner/*.yml")
FRAG_SDF_LIBRARIES = os.path.join(test_path, "frag/input_lib_sdf.yaml")
FRAG_PDB_LIBRARIES = os.path.join(test_path, "frag/input_lib_pdb.yaml")
FRAG_ANALYSIS_TO_POINT = os.path.join(test_path, "frag/input_point_analysis.yaml")
FRAGMENT_ATOM = os.path.join(test_path, "frag/input_fragment_atom.yaml")
FRAG_WATERS = os.path.join(test_path, "frag/input_frag_waters.yaml")

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

water_lines = [
     'HETATM 2538  OW  HOH B 162      51.000  92.000  14.000  1.00  0.00           O',
     'HETATM 2539 1HW  HOH B 162      51.757  92.586  14.000  1.00  0.00           H',
     'HETATM 2540 2HW  HOH B 162      50.243  92.586  14.000  1.00  0.00           H',
     'HETATM 2541  OW  HOH B 163      71.000  60.000  20.000  1.00  0.00           O',
     'HETATM 2542 1HW  HOH B 163      71.757  60.586  20.000  1.00  0.00           H',
     'HETATM 2543 2HW  HOH B 163      70.243  60.586  20.000  1.00  0.00           H',
     'HETATM 2544  OW  HOH B 164      81.000  89.000  87.000  1.00  0.00           O',
     'HETATM 2545 1HW  HOH B 164      81.757  89.586  87.000  1.00  0.00           H',
     'HETATM 2546 2HW  HOH B 164      80.243  89.586  87.000  1.00  0.00           H',
     'HETATM 2547  OW  HOH B 165      21.000  91.000  48.000  1.00  0.00           O',
     'HETATM 2548 1HW  HOH B 165      21.757  91.586  48.000  1.00  0.00           H',
     'HETATM 2549 2HW  HOH B 165      20.243  91.586  48.000  1.00  0.00           H'

]


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
        if os.path.exists(path):
            shutil.rmtree(path, ignore_errors=True)

    main.run_platform_from_yaml(ext_args)
    captured = capsys.readouterr()

    new_output_path = glob.glob(output)[0]
    top_results = glob.glob(os.path.join(new_output_path, "top_result", "*pdb"))

    assert "Skipped - FragPELE will not run." not in captured.out
    assert os.path.exists(new_output_path)
    assert len(top_results) == 2


def test_flags(capsys, ext_args=FLAGS_ARGS, output="water_processed_aminoCA1N1"):
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
    job = main.run_platform_from_yaml(ext_args)
    captured = capsys.readouterr()

    errors = tests.utils.check_file(
        job.working_dir[0],
        "control_folder/0_pele_template.conf",
        td.PELE_VALUES + FRAG_FLAGS,
        errors,
    )
    errors = tests.utils.check_file(
        job.working_dir[0], "DataLocal/LigandRotamerLibs/SB4.rot.assign", "60", errors
    )


    assert not errors
    assert "Skipped - FragPELE will not run." not in captured.out



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
        main.run_platform_from_yaml(file)


@pytest.mark.parametrize(
    ("yaml_file", "expected_lines"),
    [(FRAG_SDF_LIBRARIES, SDF_lines), (FRAG_PDB_LIBRARIES, PDB_lines)],
)
def test_libraries(capsys, yaml_file, expected_lines):
    """
    Tests the growing of fragments from a custom-made SDF and PDB libraries.
    Tests the asymmetric hydrogen detector.

    Parameters
    ----------
    yaml_file : str
        Path to PELE input file.
    expected_lines : list[str]
        List of lines expected in input.conf.
    """
    errors = []
    if os.path.exists("input.conf"):
        os.remove("input.conf")

    job = main.run_platform_from_yaml(yaml_file)
    captured = capsys.readouterr()

    for path in job.working_dir:
        analysis_folder = os.path.join(path, "results")
        assert os.path.isdir(analysis_folder)

    assert "Skipped - FragPELE will not run." not in captured.out
    assert not errors


def test_analysis_to_point(ext_args=FRAG_ANALYSIS_TO_POINT):
    """
    Tests the automated analysis to retrieve most promising fragments
    from a custom-made library based on their proximity to a certain point.
    Tests ligand clustering.

    Parameters
    ----------
    ext_args : str
        Path to PELE input file.
    """
    main.run_platform_from_yaml(ext_args)

    errors = tests.utils.check_file(
        os.getcwd(), "point_analysis.csv", point_analysis_lines, [], ",", 4
    )

    assert not errors


def test_fragment_atom(capsys, ext_args=FRAGMENT_ATOM):
    """
    Tests the frag_core_atom flag.
    
    Parameters
    ----------
    ext_args : str
        Path to PELE input file.
    """
    if os.path.exists("input.conf"):
        os.remove("input.conf")

    main.run_platform_from_yaml(ext_args)
    captured = capsys.readouterr()
    assert "Skipped - FragPELE will not run." not in captured.out


def test_frag_waters(ext_args=FRAG_WATERS):
    """
    Check if water molecules are added to the system.
  
    Parameters
    ----------
    ext_args : str
        Path to PELE input file,
    """
    water_output = []
    job = main.run_platform_from_yaml(ext_args)
    for path in job.working_dir:
        output = glob.glob(os.path.join(path, "*_top.pdb"))[0]

        with open(output, "r") as file:
            for line in file.readlines():
                if line[17:21].strip() == "HOH":
                    water_output.append(line.strip())

    assert water_lines == water_output
