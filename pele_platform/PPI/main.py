import os
import glob
from pele_platform.site_finder.cluster import cluster_best_structures
from pele_platform.PPI.simulation_launcher import launch_simulation
from pele_platform.PPI.preparation import prepare_structure
from pele_platform.Utilities.Helpers.helpers import cd, is_repeated, is_last
import pele_platform.Utilities.Parameters.parameters as pv
import pele_platform.Adaptive.simulation as si


def run_ppi(parsed_yaml: dict) -> (pv.ParametersBuilder, pv.ParametersBuilder):


    # Let user choose working folder
    original_dir = os.path.abspath(os.getcwd())
    working_folder = os.path.abspath("{}_Pele".format(parsed_yaml.residue))
    if not parsed_yaml.folder:
        working_folder = is_repeated(working_folder) if not parsed_yaml.adaptive_restart else is_last(working_folder)
    else:
        working_folder = os.path.abspath(parsed_yaml.folder)

    # Set main folder
    parsed_yaml.folder = os.path.join(working_folder, "1_interface_exploration")

    # Check n_waters before launching the simulation
    #water_checker(parsed_yaml)

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
    simulation1 = si.run_adaptive(parsed_yaml)
    simulation1_path = os.path.join(simulation1.pele_dir, simulation1.output)

    # cluster best structures
    if not parsed_yaml.debug:
        with cd(simulation1_path):
            cluster_best_structures("5", n_components=simulation1.n_components,
                residue=simulation1.residue, topology=simulation1.topology,
                directory=working_folder, logger=simulation1.logger)
    
    # adjust original input.yaml
    if not parsed_yaml.skip_refinement:
        parsed_yaml.system = os.path.join(working_folder, "refinement_input/*.pdb")
        parsed_yaml.folder = os.path.join(working_folder, "2_refinement_simulation")
        parsed_yaml.induced_fit_exhaustive = None
        parsed_yaml.ppi = None
        parsed_yaml.poses = None
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
        with cd(original_dir):
            if not parsed_yaml.debug:
                simulation2 = launch_simulation(parsed_yaml)
            else:
                simulation2 = None
    else:
        simulation2 = None
    return simulation1, simulation2
