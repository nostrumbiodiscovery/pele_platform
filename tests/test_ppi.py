from pele_platform.PPI.main import run_ppi
from pele_platform.PPI.preparation import prepare_structure 
from pele_platform.main import parseargs_yaml
from pele_platform.Utilities.Helpers.yaml_parser import YamlParser
from pele_platform.constants import constants as cs
from pele_platform import main
import pandas as pd
import glob
import os


test_path = os.path.join(cs.DIR, "Examples")
yaml = os.path.join(test_path, "PPI/input.yaml")


def test_ppi(energy_result=-2.18, yaml=yaml):
  
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


def test_prepare_structure():

    protein_file = os.path.join(test_path, "PPI/1tnf_prep.pdb")
    new_protein_file = "1tnf_prep_prep.pdb"
    ligand_pdb = os.path.join(test_path, "PPI/1tnf_ligand.pdb")
    chain = ["A", "B"]

    prepare_structure(protein_file, ligand_pdb, chain)

    with open(new_protein_file, "r") as file:
        lines = file.readlines()
        chains = []
        for line in lines:
            if line.startswith("HETATM") or line.startswith("ATOM"):
                chains.append(line[21:22].strip())
    
    assert "C" not in chains
    
    for f in glob.glob("*_prep.pdb"):
        os.remove(f)
