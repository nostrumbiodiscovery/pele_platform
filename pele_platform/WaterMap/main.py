import os

# PELE imports
import pele_platform.Adaptive.simulation as sim
import pele_platform.WaterMap.preparation as prep
import pele_platform.WaterMap.analysis as an


def run_watermap(parsed_yaml):

    # Remove all waters from the system
    user_radius = parsed_yaml.water_radius if parsed_yaml.water_radius else 6.0
    parsed_yaml.system = prep.prepare_system(parsed_yaml.system, parsed_yaml.water_center, user_radius)

    # Launch adaptive simulation
    simulation = sim.run_adaptive(parsed_yaml)

    # Get path to simulation output
    simulation_output = os.path.join(simulation.pele_dir, simulation.output)
    
    # Analyse
    analysis = an.main(simulation.water_center, simulation.water_radius, simulation_output)

    return analysis
