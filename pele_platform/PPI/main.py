from pele_platform.PPI.cluster import cluster_best_structures
from pele_platform.PPI.simulation_launcher import launch_global_exploration, launch_refinement
import yaml
from pele_platform.Utilities.Helpers.helpers import cd
import os

def run_ppi(parsed_yaml):

    # start initial simulation

    parsed_yaml.full = True
    simulation = launch_global_exploration(parsed_yaml)
    simulation_path = os.path.join(simulation.pele_dir, simulation.output)
    
    # get best structures and cluster them
    with cd(simulation_path):
        # NEED ALGORITHM TO CHOOSE OPTIMUM NUMBERS OF CLUSTERS!!!!
        #cluster_best_structures("5", n_components = parsed_yaml.cpus-1)
        cluster_best_structures("5", n_components=simulation.n_components,
            residue=simulation.residue, topology=simulation.topology)
    
    # adjust original input.yaml
    parsed_yaml.system = os.path.join(simulation_path, "refinement_input/*.pdb")
    parsed_yaml.folder = "refinement_simulation"
    parsed_yaml.full = None
    parsed_yaml.poses = None
    parsed_yaml.induced_fit_exhaustive = True
    parsed_yaml.box_center = simulation.box_center
    parsed_yaml.box_radius = simulation.box_radius
        
    # refine selected best structures
    with cd(simulation.pele_dir):
    	induced_fit = launch_refinement(parsed_yaml)

    return simulation, induced_fit
