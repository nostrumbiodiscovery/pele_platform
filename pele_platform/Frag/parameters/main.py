from pele_platform.Frag.parameters import  files, simulation, opcionals, water, metrics
import pele_platform.Utilities.Parameters.pele_env as pele



class FragParameters(pele.EnviroBuilder, water.FragWaterParams,
    files.FragInputFiles, simulation.FragSimulationParameters,
    metrics.FragMetrics):

    def __init__(self, args):
        self.software = "Frag"

        # Platform common variables
        self.build_frag_variables(args)

        #Water Parameters
        water.FragWaterParams.__init__(self, args.frag_core)

        #Frag Input_files parameters
        files.FragInputFiles.__init__(self, args)
        
        # Frag Simulation Parameters
        simulation.FragSimulationParameters.__init__(self, args)

        #Frag Metric Parameters
        metrics.FragMetrics.__init__(self, args)

        # Frag Optional Parameters
        opcionals.FragOpcionalParameters.__init__(self, args)

        # Keep initial arguments
        self.args = args



