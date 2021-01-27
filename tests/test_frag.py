import os
import glob
import pytest
import shutil
import pele_platform.constants.constants as cs
import pele_platform.main as main
import pele_platform.Frag.checker as ch
import pele_platform.Errors.custom_errors as ce
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

PDB_lines = [
    "../pele_platform/Examples/frag/pdb_test_lib/mol1.pdb C3-H2 C1-H1",
    "../pele_platform/Examples/frag/pdb_test_lib/mol1.pdb C3-H2 C1-H2",
    "../pele_platform/Examples/frag/pdb_test_lib/mol1.pdb C3-H2 C1-H3",
    "../pele_platform/Examples/frag/pdb_test_lib/mol1.pdb C3-H2 C2-H4",
    "../pele_platform/Examples/frag/pdb_test_lib/mol1.pdb C3-H2 C2-H5",
    "../pele_platform/Examples/frag/pdb_test_lib/mol1.pdb C3-H2 N1-H6",
    "../pele_platform/Examples/frag/pdb_test_lib/mol1.pdb C3-H2 N1-H7",
    "../pele_platform/Examples/frag/pdb_test_lib/mol2.pdb C3-H2 C1-H1",
    "../pele_platform/Examples/frag/pdb_test_lib/mol2.pdb C3-H2 C1-H2",
    "../pele_platform/Examples/frag/pdb_test_lib/mol2.pdb C3-H2 C1-H3",
    "../pele_platform/Examples/frag/pdb_test_lib/mol2.pdb C3-H2 N1-H4",
    "../pele_platform/Examples/frag/pdb_test_lib/mol2.pdb C3-H2 N1-H5",
]

SDF_lines = [
    "../pele_platform/Examples/frag/sdf_test_lib/test_library-1.pdb C3-H2 C1-H1",
    "../pele_platform/Examples/frag/sdf_test_lib/test_library-1.pdb C3-H2 C1-H2",
    "../pele_platform/Examples/frag/sdf_test_lib/test_library-1.pdb C3-H2 C1-H3",
    "../pele_platform/Examples/frag/sdf_test_lib/test_library-1.pdb C3-H2 C2-H4",
    "../pele_platform/Examples/frag/sdf_test_lib/test_library-1.pdb C3-H2 C2-H5",
    "../pele_platform/Examples/frag/sdf_test_lib/test_library-1.pdb C3-H2 N1-H6",
    "../pele_platform/Examples/frag/sdf_test_lib/test_library-1.pdb C3-H2 N1-H7",
    "../pele_platform/Examples/frag/sdf_test_lib/test_library-2.pdb C3-H2 C1-H1",
    "../pele_platform/Examples/frag/sdf_test_lib/test_library-2.pdb C3-H2 C1-H2",
    "../pele_platform/Examples/frag/sdf_test_lib/test_library-2.pdb C3-H2 C1-H3",
    "../pele_platform/Examples/frag/sdf_test_lib/test_library-2.pdb C3-H2 N1-H4",
    "../pele_platform/Examples/frag/sdf_test_lib/test_library-2.pdb C3-H2 N1-H5",
]

point_analysis_lines = [
    "../pele_platform/Examples/frag/analysis_data/1w7h_preparation_structure_2w_processed_mol1C3-H2C1-H1,../pele_platform/Examples/frag/analysis_data/1w7h_preparation_structure_2w_processed_mol1C3-H2C1-H1/top_result/epochsampling_result_trajectory_1.1_BindingEnergy-23.4636.pdb,2.73029273852075,../pele_platform/Examples/frag/analysis_data/1w7h_preparation_structure_2w_processed_mol1C3-H2C1-H1/top_result/epochsampling_result_trajectory_2.1_BindingEnergy-25.1634.pdb,-25.1634,../pele_platform/Examples/frag/analysis_data/1w7h_preparation_structure_2w_processed_mol1C3-H2C1-H1/top_result/epochsampling_result_trajectory_1.1_BindingEnergy-23.4636.pdb,0.7087330424726306,2.730292738520753,-23.4636"
]


@pytest.mark.parametrize("yaml", [FRAG_SIM_ARGS, FRAG_CORE_ARGS])
def test_frag_simulation(yaml):
    """
    Run end to end FragPELE test simulation starting both from input.conf and an SD file with fully grown ligands.
    """
    output = "1w7h_preparation_structure_2w_aminoC1N1"

    if os.path.exists(output):
        shutil.rmtree(output)
    job = main.run_platform(yaml)


def test_flags(ext_args=FLAGS_ARGS, output="water_processed_aminoCA1N1"):
    FRAG_FLAGS = [
        '"seed" : 3000',
    ]
    errors = []
    if os.path.exists(output):
        shutil.rmtree(output, ignore_errors=True)
    job = main.run_platform(ext_args)
    folder = output
    # if not os.path.exists(os.path.join(folder, "DataLocal/LigandRotamerLibs/STR.rot.assign")) or not os.path.exists(os.path.join(folder, "DataLocal/LigandRotamerLibs/MG.rot.assign")):
    # errors.append("External rotamer flag not working")
    # if not os.path.exists(os.path.join(folder, "DataLocal/Templates/OPLS2005/HeteroAtoms/strz")) or not os.path.exists(os.path.join(folder, "DataLocal/Templates/OPLS2005/HeteroAtoms/mgz")):
    # errors.append("External templates flag not working")
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
    files = glob.glob(ext_args)
    for file in files:
        try:
            job = main.run_platform(file)
        except Exception:
            assert False


@pytest.mark.parametrize(
    ("yaml", "expected"),
    [(FRAG_SDF_LIBRARIES, SDF_lines), (FRAG_PDB_LIBRARIES, PDB_lines)],
)
def test_libraries(yaml, expected):
    """
    Runs FragPELE in debug mode to check if the input.conf file is correctly enumerated based on available SDF or PDB
    fragment library.
    """
    if os.path.exists("input.conf"):
        os.remove("input.conf")

    job = main.run_platform(yaml)
    errors = []
    errors = td.check_file(os.getcwd(), "input.conf", expected, errors)
    assert not errors

    # remove PDB files created during conversion from SDF
    temp_files = glob.glob(os.path.join(test_path, "frag/sdf_test_lib/*.pdb"))
    for f in temp_files:
        if f:
            os.remove(f)


def test_analysis_to_point(ext_args=FRAG_ANALYSIS_TO_POINT):
    """
    Runs FragPELE in only_analysis mode to check, if the software correctly analysis mock simulation output.
    """
    job = main.run_platform(ext_args)
    errors = []
    errors = td.check_file(
        os.getcwd(), "point_analysis.csv", point_analysis_lines, errors
    )
    assert not errors


def test_checker_nlimits():
    sdf = os.path.join(test_path, "frag/checker/receptor.sdf")
    try:
        ch.check_limit_number_atoms(sdf, 100)
    except ce.LigandSizeExceed:
        assert True
    else:
        assert False


def test_checker_subsearch():
    from rdkit import Chem

    ligand = os.path.join(test_path, "frag/checker/4RFM_series.sdf")
    core = os.path.join(test_path, "frag/checker/4RFM_proc.pdb")

    mol = next(Chem.SDMolSupplier(ligand))
    core = Chem.rdmolfiles.MolFromPDBFile(core)
    atoms_in_common = mol.GetSubstructMatches(core)[0]
    atoms_in_common_after = ch.check_substructure_match(core, mol, atoms_in_common)

    assert atoms_in_common != atoms_in_common_after
    assert atoms_in_common_after[atoms_in_common.index(13)] == 12
    assert atoms_in_common_after[atoms_in_common.index(12)] == 13
