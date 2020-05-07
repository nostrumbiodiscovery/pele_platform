from pele_platform.Frag.parameters import  files, simulation, opcionals
import pele_platform.Utilities.Parameters.pele_env as pele


class FragParameters(pele.EnviroBuilder, 
    files.FragInputFiles, simulation.FragSimulationParameters):

    def __init__(self, args):
        self.software = "Frag"

        #Platform common variables
        self.build_frag_variables(args)

        #Frag Input_files parameters
        files.FragInputFiles.__init__(self, args)
        
        #Frag Simulation Parameters
        simulation.FragSimulationParameters.__init__(self, args)

        #Frag Opcional Parameters
        opcionals.FragOpcionalParameters.__init__(self, args)

        #Keep inital arguments
        self.args = args



