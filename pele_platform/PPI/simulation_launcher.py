import os
import pele_platform.Adaptive.simulation as si

DIR = os.path.dirname(os.path.abspath(__file__))


def launch_simulation(original_yaml):
    job_parameters = si.run_adaptive(original_yaml)
    return job_parameters
