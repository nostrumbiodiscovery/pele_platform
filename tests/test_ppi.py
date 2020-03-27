from pele_platform.PPI.main import run_ppi
from pele_platform.main import parseargs_yaml, YamlParser
import pandas as pd
import glob
import os

def test_PPI(energy_result=-3.87):

    # parse input.yaml
    original_yaml = os.path.join(os.getcwd(), "data/input_global.yaml")
    arguments = parseargs_yaml([original_yaml,])
    arguments = YamlParser(arguments.input_file)

    # run the platform
    run_ppi(arguments)

    # checkpoints
    output_csv = pd.read_csv("global_simulation/output/clustering_output.csv")
    best_energy = round(output_csv["binding_energy"].min(),2)
    nfiles = len(glob.glob("global_simulation/output/refinement_input/*.pdb"))
    nfiles_refinement = len(glob.glob("global_simulation/refinement_simulation/results/BestStructs/epoch*"))

    # test
    assert nfiles == arguments.n_components 
    assert best_energy == energy_result
    assert nfiles_refinement
