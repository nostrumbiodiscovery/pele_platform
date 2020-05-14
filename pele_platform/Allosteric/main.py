from pele_platform.Allosteric.cluster import cluster_best_structures
from pele_platform.Allosteric.simulation_launcher import launch_global_exploration, launch_refinement
from pele_platform.Utilities.Helpers.helpers import cd, is_repited, is_last
import os


def run_allosteric(parsed_yaml):
    #Let user choose working folder
    original_dir = os.path.abspath(os.getcwd())
    working_folder = os.path.abspath("{}_Pele".format(parsed_yaml.residue))
    if not parsed_yaml.folder:
        working_folder = is_repited(working_folder) if not parsed_yaml.adaptive_restart else is_last(working_folder)
    else:
        working_folder = os.path.abspath(parsed_yaml.folder)

    #Set main folder
    parsed_yaml.folder = os.path.join(working_folder, "1_global_exploration")

    # start initial simulation
    parsed_yaml.full = True
    simulation = launch_global_exploration(parsed_yaml)
    simulation_path = os.path.join(simulation.pele_dir, simulation.output)
    
    # get best structures and cluster them
    with cd(simulation_path):
        # NEED ALGORITHM TO CHOOSE OPTIMUM NUMBERS OF CLUSTERS!!!!
        cluster_best_structures("5", n_components=simulation.n_components,
            residue=simulation.residue, topology=simulation.topology,
            directory=working_folder)
    
    # adjust original input.yaml
    if not parsed_yaml.skip_refinement:
        parsed_yaml.system = os.path.join(working_folder, "refinement_input/*.pdb")
        parsed_yaml.folder = os.path.join(working_folder, "2_refinement_simulation")
        parsed_yaml.full = None
        parsed_yaml.poses = None
        parsed_yaml.induced_fit_exhaustive = True
        if not parsed_yaml.test:
            parsed_yaml.iterations = 100
            parsed_yaml.pele_steps = 10
        parsed_yaml.box_center = simulation.box_center
        parsed_yaml.box_radius = simulation.box_radius
        # refine selected best structures
        with cd(original_dir):
            induced_fit = launch_refinement(parsed_yaml)
    else:
        induced_fit = None

    return simulation, induced_fit
