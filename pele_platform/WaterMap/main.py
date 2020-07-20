import os

# PELE imports
import pele_platform.Adaptive.simulation as sim
import pele_platform.WaterMap.preparation as prep
import pele_platform.WaterMap.analysis as an


def run_watermap(parsed_yaml):

    # Remove all waters from the system
    parsed_yaml.system = prep.remove_water(parsed_yaml.system)

    # Launch adaptive simulation
    simulation = sim.run_adaptive(parsed_yaml)

    # Get path to simulation output
    simulation_output = os.path.join(simulation.pele_dir, simulation.output)
    
    # Analyse
    analysis = an.Watermap(simulation_output)
    analysis_results = analysis.run()

    return analysis_results
