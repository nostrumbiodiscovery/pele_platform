from pele_platform.Allosteric.cluster import cluster_best_structures
from pele_platform.PPI.simulation_launcher import launch_induced
from pele_platform.PPI.preparation import prepare_structure
import yaml
from pele_platform.Utilities.Helpers.helpers import cd
import os

def run_ppi(parsed_yaml):

    # get arguments from input.yaml
    protein_file = parsed_yaml.system
    chain = parsed_yaml.protein
    ligand_pdb = parsed_yaml.ligand_pdb
   
    # remove chains except for "protein" flag
    protein_file = prepare_structure(protein_file, ligand_pdb, chain)

    # start simualtion 1 - induced fit
    parsed_yaml.induced_fit_exhaustive = True
    parsed_yaml.randomize = True
    simulation1 = launch_induced(parsed_yaml)
    simulation1_path = os.path.join(simulation1.pele_dir, simulation1.output)
    
    # cluster best structures
    with cd(simulation1_path):
        cluster_best_structures("5", n_components=simulation1.n_components,
            residue=simulation1.residue, topology=simulation1.topology)
    
    # adjust original input.yaml
    #parsed_yaml.system = os.path.join(simulation_path, "refinement_input/*.pdb")
    #parsed_yaml.folder = "refinement_simulation"
    #parsed_yaml.full = None
    #parsed_yaml.poses = None
    #parsed_yaml.induced_fit_exhaustive = True
    #parsed_yaml.box_center = simulation.box_center
    #parsed_yaml.box_radius = simulation.box_radius
        
    # refine selected best structures
    #with cd(simulation.pele_dir):
    #	induced_fit = launch_refinement(parsed_yaml)

    #return simulation
