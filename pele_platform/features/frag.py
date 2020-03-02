import os
import pele_platform.constants.constants as cs
import pele_platform.constants.pele_params as pcs

"""
Description of the file: This file specifies what files to use
for each job.

"""
    
def retrieve_software_settings(args):

        SOFTWARE_CONSTANTS = {
                 "simulation_params" : {
                 "standard": {"params": pcs.FRAG, "COMligandConstraint": 1}
        }
        }
        
        software_setings = SOFTWARE_CONSTANTS
        type_simulation = "standard"
        software_setings["simulation_params"] = software_setings["simulation_params"].get(type_simulation, {})
        return software_setings
