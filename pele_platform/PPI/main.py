from pele_platform.Allosteric.cluster import cluster_best_structures
from pele_platform.PPI.simulation_launcher import launch_simulation
from pele_platform.PPI.preparation import prepare_structure
from pele_platform.Utilities.Helpers.helpers import cd
import os


def run_ppi(parsed_yaml):

    # get arguments from input.yaml
    protein_file = parsed_yaml.system
    chain = parsed_yaml.protein
    ligand_pdb = parsed_yaml.ligand_pdb
   
    # remove chains except for "protein" flag
    protein_file = prepare_structure(protein_file, ligand_pdb, chain)
    parsed_yaml.system = protein_file

    # start simualtion 1 - induced fit
    parsed_yaml.induced_fit_exhaustive = True
    simulation1 = launch_simulation(parsed_yaml)
    simulation1_path = os.path.join(simulation1.pele_dir, simulation1.output)
    
    # cluster best structures
    with cd(simulation1_path):
        cluster_best_structures("5", n_components=simulation1.n_components,
            residue=simulation1.residue, topology=simulation1.topology)
    
    # adjust original input.yaml
    parsed_yaml.system = os.path.join(simulation1_path, "refinement_input/*.pdb")
    parsed_yaml.folder = "refinement_simulation"
    parsed_yaml.induced_fit_exhaustive = None
    parsed_yaml.ppi = None
    parsed_yaml.poses = None
    parsed_yaml.rescoring = True
    parsed_yaml.iterations = 1
    parsed_yaml.steps = 100
    parsed_yaml.box_center = simulation1.box_center
    parsed_yaml.box_radius = 100  # We should have a look at how to set no box but at the moment super big
        
    # start simulation 2 - minimisation
    with cd(simulation1.pele_dir):
        simulation2 = launch_simulation(parsed_yaml)

    return simulation1, simulation2
