import matplotlib
matplotlib.use("Agg")
import sys
import pele_platform.constants.constants as cs
sys.path.append(cs.DIR)
from argparse import HelpFormatter
from operator import attrgetter
import argparse
import os
import pele_platform.Adaptive.simulation as ad
import pele_platform.Utilities.Helpers.yaml_parser as yp
import pele_platform.Frag.simulation as fr
import pele_platform.Utilities.Helpers.environment_variables as ev
from pele_platform.PPI.main import run_ppi

class Launcher():


    def __init__(self, arguments):
        self.cpus = arguments.cpus
        self.restart = arguments.restart
        self.test = arguments.test
        self._args = arguments
        if arguments.frag_core:
            self._args.pele_feature = "frag"
            self._args.system = self._args.frag_core
        elif arguments.ppi:
            self._args.pele_feature = "PPI"
        else: 
            self._args.pele_feature = "adaptive"

    def launch(self):
        if not self._args.no_check:
            ev.check_variables(self._args)
        if self._args.pele_feature == "adaptive":
            job_variables = ad.run_adaptive(self._args)
        elif self._args.pele_feature == "PPI":
            job_variables = run_ppi(self._args)
        elif self._args.pele_feature == "frag":
            #Set variables and input ready 
            job_variables = fr.FragRunner(self._args)
            job_variables.prepare_control_file()
            #Set test variables if desired
            if self.test:
                job_variables.set_test_variables()
            #Depending on input different method
            if job_variables.ligands: #Full ligands as sdf
                job_variables.prepare_input_file()
                job_variables.run()
            else:
                job_variables.run()
            # Execute job
        return job_variables


def parseargs_yaml(args=[]):
    parser = argparse.ArgumentParser(description='Run PELE Platform')
    parser.add_argument('input_file', type=str, help='Yaml input file')
    args = parser.parse_args(args) if args else parser.parse_args()
    return args
    
def run_platform(input_yaml):
    arguments = parseargs_yaml([input_yaml,])
    arguments = yp.YamlParser(arguments.input_file)
    job_params = Launcher(arguments).launch()
    return job_params


if __name__ == "__main__":
    arguments = parseargs_yaml()
    arguments = yp.YamlParser(arguments.input_file)
    job_params = Launcher(arguments).launch()
