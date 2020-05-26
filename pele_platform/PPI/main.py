from pele_platform.Allosteric.cluster import cluster_best_structures
from pele_platform.PPI.simulation_launcher import launch_simulation
from pele_platform.PPI.preparation import prepare_structure
from pele_platform.Utilities.Helpers.water import add_water, water_checker
from pele_platform.Utilities.Helpers.helpers import cd
import glob
import os


def run_ppi(parsed_yaml):

    # Check n_waters before launching the simulation
    water_checker(parsed_yaml)

    # get arguments from input.yaml
    n_waters = parsed_yaml.n_waters
    parsed_yaml.n_waters = None
    protein_file = parsed_yaml.system
    chain = parsed_yaml.protein
    ligand_pdb = parsed_yaml.ligand_pdb

    # no waters in the first simulation
    parsed_yaml.water_arg = None 
   
    # remove chains except for "protein" flag
    protein_file = prepare_structure(protein_file, ligand_pdb, chain, True)
    parsed_yaml.system = protein_file

    # start simulation 1 - induced fit
    parsed_yaml.induced_fit_exhaustive = True
    simulation1 = launch_simulation(parsed_yaml)
    simulation1_path = os.path.join(simulation1.pele_dir, simulation1.output)
    
    # cluster best structures
    with cd(simulation1_path):
        cluster_best_structures("5", n_components=simulation1.n_components,
            residue=simulation1.residue, topology=simulation1.topology)
   
    if not parsed_yaml.skip_refinement:

        # adjust original input.yaml
        parsed_yaml.system = os.path.join(simulation1_path, "refinement_input/*.pdb")
        parsed_yaml.folder = "refinement_simulation"
        parsed_yaml.induced_fit_exhaustive = None
        parsed_yaml.poses = None
        parsed_yaml.ppi = False
        parsed_yaml.rescoring = True
        del parsed_yaml.water_arg
        # Set waters ony if specified by user
        if n_waters != 0:
            parsed_yaml.waters = "all_waters"
            parsed_yaml.n_waters = n_waters
        else:
            parsed_yaml.waters = None
            parsed_yaml.n_waters = n_waters
        parsed_yaml.adaptive_restart = False
        if not parsed_yaml.test:
            parsed_yaml.iterations = 1
            parsed_yaml.steps = 100
    
        parsed_yaml.box_center = simulation1.box_center
        parsed_yaml.box_radius = 100  # We should have a look at how to set no box but at the moment super big
    
        # start simulation 2 - minimisation
        with cd(simulation1.pele_dir):
            simulation2 = launch_simulation(parsed_yaml)
    else:
        simulation2 = None

    return simulation1, simulation2

