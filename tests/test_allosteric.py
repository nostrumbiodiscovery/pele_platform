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

yaml = os.path.join(test_path, "Allosteric/input_skipref.yaml")
def test_allosteric_skipref(yaml=yaml):

    job, = main.run_platform(yaml)
    refinement_simulation = os.path.join(os.path.dirname(job.pele_dir), "2_Selection")
    
    assert not os.path.exists(refinement_simulation)


def test_allosteric_pdb():
    yaml = os.path.join(test_path, "Allosteric/input_global.yaml")
    job, _, _ = main.run_platform(yaml)

    # checkpoints
    refinement_input = glob.glob(os.path.join(os.path.dirname(job.pele_dir), "2_Selection/*.pdb"))
    nfiles_refinement = len(glob.glob(os.path.join(job.pele_dir, "results/BestStructs/epoch*")))
    best_energy_input = os.path.join(os.path.dirname(job.pele_dir), "2_Selection",  "epoch0_trajectory_3.1_BindingEnergy-1.5142.pdb")
    
    # test
    assert best_energy_input in refinement_input
    assert nfiles_refinement > 0


yaml = os.path.join(test_path, "Allosteric/input_global_xtc.yaml")
def test_allosteric_xtc(yaml=yaml):

    job, _, _ = main.run_platform(yaml)

    # checkpoints
    best_energy_input = os.path.join(os.path.dirname(job.pele_dir), "2_Selection",  "epoch0_trajectory_3.1_BindingEnergy-1.5142.pdb")
    refinement_input = glob.glob(os.path.join(os.path.dirname(job.pele_dir), "2_Selection/*.pdb"))
    nfiles_refinement = len(glob.glob(os.path.join(job.pele_dir, "results/BestStructs/epoch*")))

    # test
    assert best_energy_input in refinement_input 
    assert nfiles_refinement > 0


yaml = os.path.join(test_path, "Allosteric/input_restart.yaml")

@pytest.mark.xfail(machine == "10.10.2.3", reason = "Fails on nbdcalc01 due to mpi. Mock simulation folder was prepared with srun.")
def test_allosteric_restart(yaml=yaml):

    if os.path.exists('allosteric'):
        shutil.rmtree('allosteric')

    os.system("cp -r ../pele_platform/Examples/Allosteric/allosteric .")
    job1, sel, job2 = main.run_platform(yaml)

    # checkpoints
    nfiles_refinement = len(glob.glob(os.path.join("allosteric/3_InducedFitExhaustive/results/BestStructs/epoch*")))
    
    # test
    assert nfiles_refinement > 0

