import os
import pele_platform.constants.constants as cs

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

def check_variables(args):
    """
    Check all external requirements are ther
    before starting the simulation 
    """
    #Create class objects
    variables = [
    EnvVariable("pele_data", args.pele_data, os.path.join(cs.PELE, "Data"), "--pele_data /path/to/data/folder/", "export PELE=/path/to/PELE-1.X/"),
    EnvVariable("pele_documents", args.pele_documents, os.path.join(cs.PELE, "Documents"), "--pele_documents /path/to/documents/folder", "export PELE=/path/to/PELE-1.X/"),
    EnvVariable("pele_exec", args.pele_exec, os.path.join(cs.PELE, "bin/Pele_mpi"), "--pele_exec /path/to/PELE_exec", "export PELE=/path/to/PELE-1.X/"),
    EnvVariable("pele_license", args.pele_license, os.path.join(cs.PELE, "licenses"), "--pele_license /path/to/licenses", "export PELE=/path/to/PELE-1.X/"),
    EnvVariable("schrodinger", args.schrodinger, cs.SCHRODINGER, "--schrodinger /path/to/schrodinger-20XX/", "export SCHRODINGER=/path/to/schrodinger-20XX/")
    ]
    #Check validity of the objects
    for variable in variables:
        variable.is_valid()
