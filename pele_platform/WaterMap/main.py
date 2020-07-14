import os

# PELE imports
import pele_platform.Adaptive.simulation as si
import pele_platform.WaterMap.preparation as pr
import pele_platform.WaterMap.analysis as an


def run_watermap(parsed_yaml):

    # Remove all waters from the system
    protein_file = pr.remove_water(parsed_yaml.system)
    parsed_yaml.system = protein_file

    # Add one water molecule and make sure it gets perturbed
    parsed_yaml.n_waters = 1
    parsed_yaml.all_waters = True

    # Launch adaptive simulation
    simulation = si.run_adaptive(parsed_yaml)

    # Get path to simulation output and cluster
    simulation_output = os.path.join(simulation.pele_dir, simulation.output)
    analysis_results = an.main(simulation_output)

    return analysis_results
