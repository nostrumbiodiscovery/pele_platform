import glob
import os

import pele_platform.Utilities.Helpers.helpers as hp
import pele_platform.Adaptive.simulation as sim
import pele_platform.WaterMap.preparation as prep
from pele_platform.analysis.water_analysis import WaterAnalysis


def run_watermap(parsed_yaml):
    """
    Runs the whole water map package and performs analysis.
    TODO: Rewrite this as a class.
    Parameters
    ----------
    parsed_yaml : YamlParser
        Arguments from input.yaml.
    Returns
    -------
        Output of the analysis.
    """

    # Remove all waters from the system
    user_radius = parsed_yaml.water_radius if parsed_yaml.water_radius else 6.0

    try:
        water_center = hp.get_coords_from_residue(
            parsed_yaml.system, parsed_yaml.water_center
        )
    except:  # TODO: Fix base exception!
        water_center = parsed_yaml.water_center
    parsed_yaml.system = prep.prepare_system(
        parsed_yaml.system, water_center, user_radius
    )
    parsed_yaml.water_center = water_center

    if not parsed_yaml.folder:
        parsed_yaml.folder = hp.get_next_peledir("Water_Pele")

    # Launch adaptive simulation
    simulation = sim.run_adaptive(parsed_yaml)

    # analyse the simulation
    results = os.path.join(simulation.pele_dir, "results")
    results_folders = glob.glob(results.replace("results", "results*"))
    breakpoint()
    if len(results_folders) == 1:
        water_analysis_folder = results_folders
    else:
        index = max(
            [
                int(os.path.basename(folder).split("_")[-1])
                for folder in results_folders
                if "_" in os.path.basename(folder)
            ]
        )
        water_analysis_folder = os.path.join(
            os.path.dirname(results),
            os.path.basename(results) + "_" + str(index),
        )

    analysis = WaterAnalysis.from_parameters(simulation)
    analysis.analyse_waters(os.path.join(water_analysis_folder, "water_analysis"))

    return simulation
