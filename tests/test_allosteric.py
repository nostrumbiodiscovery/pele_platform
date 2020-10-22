import shutil
import pandas as pd
import glob
import os
from pele_platform.constants import constants as cs
from pele_platform import main

test_path = os.path.join(cs.DIR, "Examples")

yaml = os.path.join(test_path, "Allosteric/input_skipref.yaml")
def test_allosteric_skipref(yaml=yaml):

    job, job2 = main.run_platform(yaml)
    refinement_simulation = os.path.join(os.path.dirname(job.pele_dir), "2_refinement_simulation")
    
    assert not os.path.exists(refinement_simulation)
    assert not job2


def test_allosteric_pdb():
    yaml = os.path.join(test_path, "Allosteric/input_global.yaml")
    job, job2 = main.run_platform(yaml)

    # checkpoints
    refinement_input = glob.glob(os.path.join(os.path.dirname(job.pele_dir), "refinement_input/*.pdb"))
    nfiles_refinement = len(glob.glob(os.path.join(job2.pele_dir, "results/BestStructs/epoch*")))
    best_energy_input = os.path.join(os.path.dirname(job.pele_dir), "refinement_input",  "epoch0_trajectory_3.1_BindingEnergy-1.5142.pdb")
    
    # test
    assert best_energy_input in refinement_input
    assert len(refinement_input) == 3
    assert nfiles_refinement > 0


yaml = os.path.join(test_path, "Allosteric/input_global_xtc.yaml")
def test_allosteric_xtc(yaml=yaml):

    job, job2 = main.run_platform(yaml)

    # checkpoints
    best_energy_input = os.path.join(os.path.dirname(job.pele_dir), "refinement_input",  "epoch0_trajectory_3.1_BindingEnergy-1.5142.pdb")
    refinement_input = glob.glob(os.path.join(os.path.dirname(job.pele_dir), "refinement_input/*.pdb"))
    nfiles_refinement = len(glob.glob(os.path.join(job2.pele_dir, "results/BestStructs/epoch*")))

    # test
    assert best_energy_input in refinement_input 
    assert len(refinement_input) == 3
    assert nfiles_refinement > 0


yaml = os.path.join(test_path, "Allosteric/input_folder.yaml")
def test_working_folder(yaml=yaml, output="allosteric_folder"):
    if os.path.exists(output): shutil.rmtree(output) 
    job, _ = main.run_platform(yaml)
    assert os.path.exists(job.folder)
    
    
#yaml = os.path.join(test_path, "Allosteric/input_restart.yaml")
#def test_allosteric_restart(yaml=yaml):
#    job, job2 = main.run_platform(yaml)
#    # checkpoints
#    output_csv = pd.read_csv(os.path.join(job.pele_dir, "output/clustering_output.csv"))
#    nfiles = len(glob.glob(os.path.join(os.path.dirname(job.pele_dir), "refinement_input/*.pdb")))
#    nfiles_refinement = len(glob.glob(os.path.join(job2.pele_dir, "results/BestStructs/epoch*")))
#    # test
#    assert nfiles == job.n_components 
#    assert nfiles_refinement

