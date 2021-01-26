import shutil
import pandas as pd
import glob
import os
import shutil
from pele_platform.constants import constants as cs
from pele_platform import main
import pytest
import socket

machine = socket.gethostbyname(socket.gethostname())
test_path = os.path.join(cs.DIR, "Examples")


def test_allosteric_skipref():
    """
    Runs allosteric package with 'skip_refinement' flag and makes sure the simulation stops after GlobalExploration
    (i.e. no 2_Selection, 3_induced... or other folders are created).
    """
    yaml = os.path.join(test_path, "Allosteric/input_skipref.yaml")
    (job,) = main.run_platform(yaml)
    refinement_simulation = os.path.join(os.path.dirname(job.pele_dir), "2_Selection")

    assert not os.path.exists(refinement_simulation)


@pytest.mark.parametrize(
    "yaml",
    [
        os.path.join(test_path, "Allosteric/input_global.yaml"),
        os.path.join(test_path, "Allosteric/input_global_xtc.yaml"),
    ],
)
def test_allosteric_pdb(yaml):
    """
    Runs end to end allosteric package with both PDB and XTC to ensure that both extensions are fully supported.
    Checks clustering output to ensure we did not break the methodology and get correct binding energies.
    Asserts the pipeline runs to the end by checking the presence of rescoring output.
    """
    output = main.run_platform(yaml)
    job = output[0]  # taking 1st out of 5 BB
    rescoring = output[-1]

    refinement_input = glob.glob(
        os.path.join(os.path.dirname(job.pele_dir), "2_Selection/*.pdb")
    )
    nfiles_refinement = len(
        glob.glob(os.path.join(job.pele_dir, "results/BestStructs/epoch*"))
    )
    best_energy_input = os.path.join(
        os.path.dirname(job.pele_dir),
        "2_Selection",
        "epoch0_trajectory_3.1_BindingEnergy-1.5142.pdb",
    )
    rescoring_output = glob.glob(
        os.path.join(rescoring.pele_dir, rescoring.output, "*/*.pdb")
    )

    assert best_energy_input in refinement_input
    assert nfiles_refinement > 0
    assert len(rescoring_output) > 1


@pytest.mark.xfail(
    machine == "10.10.2.3",
    reason="Fails on nbdcalc01 due to mpi. Mock simulation folder was prepared with srun.",
)
def test_allosteric_restart():
    yaml = os.path.join(test_path, "Allosteric/input_restart.yaml")
    if os.path.exists("allosteric"):
        shutil.rmtree("allosteric")

    os.system("cp -r ../pele_platform/Examples/Allosteric/allosteric .")
    output = main.run_platform(yaml)
    nfiles_refinement = len(
        glob.glob(
            os.path.join("allosteric/3_InducedFitExhaustive/results/BestStructs/epoch*")
        )
    )
    assert nfiles_refinement > 0
