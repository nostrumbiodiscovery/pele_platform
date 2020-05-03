from pele_platform.PPI.main import run_ppi
from pele_platform.main import parseargs_yaml
from pele_platform.Utilities.Helpers.yaml_parser import YamlParser
from pele_platform.constants import constants as cs
from pele_platform import main
import pandas as pd
import glob
import os

test_path = os.path.join(cs.DIR, "Examples")
yaml = os.path.join(test_path, "ppi/input_global.yaml")

def test_PPI(energy_result=-3.26, yaml=yaml):
    #Function to test
    job, _ = main.run_platform(yaml)

    # checkpoints
    output_csv = pd.read_csv(os.path.join(job.pele_dir, "output/clustering_output.csv"))
    best_energy = round(output_csv["binding_energy"].min(),2)
    nfiles = len(glob.glob(os.path.join(job.pele_dir, "output/refinement_input/*.pdb")))
    nfiles_refinement = len(glob.glob(os.path.join(job.pele_dir, "refinement_simulation/results/BestStructs/epoch*")))

    # test
    assert nfiles == job.n_components 
    assert best_energy == energy_result
    assert nfiles_refinement

yaml = os.path.join(test_path, "ppi/input_global_xtc.yaml")
def test_PPI_xtc(energy_result=-3.26, yaml=yaml):
    #Function to test
    job, _ = main.run_platform(yaml)

    # checkpoints
    output_csv = pd.read_csv(os.path.join(job.pele_dir, "output/clustering_output.csv"))
    best_energy = round(output_csv["binding_energy"].min(),2)
    nfiles = len(glob.glob(os.path.join(job.pele_dir, "output/refinement_input/*.pdb")))
    nfiles_refinement = len(glob.glob(os.path.join(job.pele_dir, "refinement_simulation/results/BestStructs/epoch*")))

    # test
    assert nfiles == job.n_components 
    assert best_energy == energy_result
    assert nfiles_refinement
