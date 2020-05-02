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
            self._check_variables()
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
            elif job_variables.ai:
                job_variables.grow_ai()
            else:
                job_variables.run()
            # Execute job
        return job_variables

    def _check_variables(self):
        variables = [ 
        EnvVariable("pele_data", self._args.pele_data, os.path.join(cs.PELE, "Data"), "--pele_data /path/to/data/folder/", "export PELE=/path/to/PELE-1.X/"),
        EnvVariable("pele_documents", self._args.pele_documents, os.path.join(cs.PELE, "Documents"), "--pele_documents /path/to/documents/folder", "export PELE=/path/to/PELE-1.X/"),
        EnvVariable("pele_exec", self._args.pele_exec, os.path.join(cs.PELE, "bin/Pele_mpi"), "--pele_exec /path/to/PELE_exec", "export PELE=/path/to/PELE-1.X/"),
        EnvVariable("pele_license", self._args.pele_license, os.path.join(cs.PELE, "licenses"), "--pele_license /path/to/licenses", "export PELE=/path/to/PELE-1.X/"),
        EnvVariable("schrodinger", self._args.schrodinger, cs.SCHRODINGER, "--schrodinger /path/to/schrodinger-20XX/", "export SCHRODINGER=/path/to/schrodinger-20XX/")
        ]
        for variable in variables:
            variable.is_valid() 


class EnvVariable():

    def __init__(self, name, variable, default, flag, env_var):
        self.name = name
        self.variable = variable
        self.default = default
        self.flag = flag
        self.env_var = env_var
        
    def is_valid(self):
        if self.variable:
            if os.path.exists(self.variable):
                return True
        elif self.default:
           if os.path.exists(self.default):
               return True

        raise ValueError("{} not found. If you have the standard installation \
export the environment variable by doing: {}.\
 Else define the location of the library via the flag: {}".format(self.name, self.env_var, self.flag))


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
