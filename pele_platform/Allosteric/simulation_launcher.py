import os
import yaml
import pele_platform.Adaptive.simulation as si

DIR = os.path.dirname(os.path.abspath(__file__))

def launch_global_exploration(original_yaml):
    job_parameters = si.run_adaptive(original_yaml)
    return job_parameters

def launch_refinement(dict, output_folder="."):

    # launch refinement using new input.yaml
    job_parameters = si.run_adaptive(dict)
    return job_parameters
