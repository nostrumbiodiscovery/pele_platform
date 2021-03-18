from pele_platform.Frag.parameters import (files, simulation, optional,
                                           water, metrics)
from pele_platform.Utilities.Parameters import parameters

# TODO deprecated, can be removed
class FragParameters(parameters.ParametersBuilder, water.FragWaterParams,
                     files.FragInputFiles, simulation.FragSimulationParameters,
                     metrics.FragMetrics):

    def __init__(self, args):
        self.software = "Frag"

        # Platform common variables
        self.build_frag_variables(args)

        # Water Parameters
        water.FragWaterParams.__init__(self, args.frag_core)

        # Frag Input_files parameters
        files.FragInputFiles.__init__(self, args)
        
        # Frag Simulation Parameters
        simulation.FragSimulationParameters.__init__(self, args)

        # Frag Metric Parameters
        metrics.FragMetrics.__init__(self, args)

        # Frag Optional Parameters
        optional.FragOptionalParameters.__init__(self, args)

        # Keep initial arguments
        self.args = args



