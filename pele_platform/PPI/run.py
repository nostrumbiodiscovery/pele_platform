from .cluster import cluster_best_structures
from .simulation_launcher import launch_global_exploration, launch_refinement
import yaml
from pele_platform.Utilities.Helpers.helpers import cd
import os

original_yaml = "input.yaml"

def run_ppi(original_yaml, ncomponents=10):

    # start initial simulation
    simulation = launch_global_exploration(original_yaml)
    simulation_path = os.path.join(simulation.pele_dir, simulation.output)
    
    # get best structures and cluster them
    with cd(simulation_path):
        cluster_best_structures("5", n_components=ncomponents)
    
    # load original input.yaml as dictionary
    with open(original_yaml) as original_yaml:
        dict = yaml.load(original_yaml, Loader=yaml.FullLoader)
    
    dict["system"] = os.path.join(simulation_path, "refinement_input/*.pdb")
    dict["working_folder"] = "refinement_simulation"
    del dict["global"]
    del dict["poses"]
    dict["out_in"] = True
    
    # refine selected best structures
    with cd(simulation.pele_dir):
    	launch_refinement(dict)
