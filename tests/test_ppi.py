import pandas as pd
import glob
import os
import shutil
from pele_platform.Utilities.BuildingBlocks.preparation import prepare_structure
from pele_platform.constants import constants as cs
from pele_platform import main

test_path = os.path.join(cs.DIR, "Examples")


def test_ppi_skipref():
    """
    Runs allosteric package with 'skip_refinement' flag and makes sure the simulation stops after the first simulation
    block (folder 3_Rescoring should not exist).
    """
    yaml = os.path.join(test_path, "PPI/input_skipref.yaml")
    (job,) = main.run_platform(yaml)
    files_refinement = glob.glob(
        os.path.join(job.pele_dir, "3_Rescoring/results/BestStructs/epoch*")
    )

    assert not files_refinement


def test_ppi():
    """
    Runs a full PPI simulation, ensures that the correct refinement inputs were selected and the rescoring output
    exists.
    """
    yaml = os.path.join(test_path, "PPI/input.yaml")
    energy_result = -2.18

    job, sel, job2 = main.run_platform(yaml)

    output_csv = pd.read_csv(os.path.join(job.pele_dir, "output/clustering_output.csv"))
    best_energy = round(output_csv["binding_energy"].min(), 2)
    nfiles = len(
        glob.glob(os.path.join(os.path.dirname(job.pele_dir), "2_Selection/*.pdb"))
    )
    nfiles_refinement = len(
        glob.glob(os.path.join(job2.pele_dir, "results/BestStructs/epoch*"))
    )

    assert nfiles == job.cpus - 1
    assert best_energy == energy_result
    assert nfiles_refinement


def test_prepare_structure():
    """
    The PPI package requires system preparation where we remove all redundant protein chains and only keep the ones
    specified by the user. This function takes a PDB file with three chains (A, B, C), preprocesses it and makes sure
    chain C and all water molecules are removed.
    """
    protein_file = os.path.join(test_path, "PPI/1tnf_prep.pdb")
    new_protein_file = "1tnf_prep_prep.pdb"
    ligand_pdb = os.path.join(test_path, "PPI/1tnf_ligand.pdb")
    chain = ["A", "B"]

    prepare_structure(protein_file, ligand_pdb, chain, remove_water=True)

    with open(new_protein_file, "r") as file:
        lines = file.readlines()
        chains = []
        no_water = True
        for line in lines:
            if "HOH" in line:
                no_water = False
            if line.startswith("HETATM") or line.startswith("ATOM"):
                chains.append(line[21:22].strip())

    assert "C" not in chains
    assert no_water

    for f in glob.glob("*_prep.pdb"):
        os.remove(f)


def test_working_folder():
    """
    Ensures setting working_folder flag works with the multi-step PPI workflow.
    """
    yaml = os.path.join(test_path, "PPI/input_folder.yaml")
    output = "ppi_folder"
    if os.path.exists(output):
        shutil.rmtree(output)
    (job,) = main.run_platform(yaml)
    assert os.path.exists(job.folder)
