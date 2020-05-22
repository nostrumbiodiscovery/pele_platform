import matplotlib
matplotlib.use("Agg")
import sys
import pele_platform.Checker.python_version as pv
pv.check_python_version()
import pele_platform.constants.constants as cs
sys.path.append(cs.DIR)
from argparse import HelpFormatter
from operator import attrgetter
import argparse
import os
import pele_platform.Adaptive.simulation as ad
import pele_platform.Utilities.Helpers.yaml_parser as yp
import pele_platform.Frag.simulation as fr
from pele_platform.Allosteric.main import run_allosteric
import pele_platform.Checker.main as ck
import pele_platform.gpcr.main as gpcr
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
        elif arguments.allosteric:
            self._args.pele_feature = "allosteric"
        elif arguments.gpcr_orth:
            self._args.pele_feature = "gpcr_orth"
        else: 
            self._args.pele_feature = "adaptive"

    def launch(self):
        if not self._args.no_check:
            ck.check_executable_and_env_variables(self._args)
        if self._args.pele_feature == "adaptive":
            job_variables = ad.run_adaptive(self._args)
        elif self._args.pele_feature == "gpcr_orth":
            job_variables = gpcr.GpcrLauncher(self._args).run_gpcr_simulation()
        elif self._args.pele_feature == "allosteric":
            job_variables = run_allosteric(self._args)
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
