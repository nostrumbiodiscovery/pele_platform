from dataclasses import dataclass
import os
import pele_platform.constants.constants as cs
from pele_platform.Checker import executable as ex
from pele_platform.Checker import environment_variables as en
import pele_platform.Utilities.Helpers.yaml_parser as yp


EXECUTABLES_TO_CHECK = ["mpirun",]


@dataclass
class Checker():

    def check_variables(self, args: yp.YamlParser) -> None:
        for env_variable in self._generate_env_variables(args):
            env_variable.is_valid()
        for executable in self._generate_executables():
            executable.is_in_path()

    def _generate_env_variables(self, args: yp.YamlParser) -> list:
        self.env_variables = [
        en.EnvVariable("pele_data", args.pele_data, os.path.join(cs.PELE, "Data"), "--pele_data /path/to/data/folder/", "export PELE=/path/to/PELE-1.X/"),
        en.EnvVariable("pele_documents", args.pele_documents, os.path.join(cs.PELE, "Documents"), "--pele_documents /path/to/documents/folder", "export PELE=/path/to/PELE-1.X/"),
        en.EnvVariable("pele_exec", args.pele_exec, os.path.join(cs.PELE, "bin/Pele_mpi"), "--pele_exec /path/to/PELE_exec", "export PELE=/path/to/PELE-1.X/"),
        en.EnvVariable("pele_license", args.pele_license, os.path.join(cs.PELE, "licenses"), "--pele_license /path/to/licenses", "export PELE=/path/to/PELE-1.X/"),
        en.EnvVariable("schrodinger", args.schrodinger, cs.SCHRODINGER, "--schrodinger /path/to/schrodinger-20XX/", "export SCHRODINGER=/path/to/schrodinger-20XX/")
        ]
        return self.env_variables


    def _generate_executables(self) -> list:
        self.executables = [ex.Executable(executable) for executable in EXECUTABLES_TO_CHECK]
        return self.executables


    
def check_executable_and_env_variables(args: yp.YamlParser):
    """
    Check all external requirements are there
    before starting the simulation.
    1) Check env variables
    2) Check executables 
    """
    checker = Checker()
    checker.check_variables(args)
