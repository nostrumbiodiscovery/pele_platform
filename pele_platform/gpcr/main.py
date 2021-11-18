from dataclasses import dataclass
import pele_platform.Adaptive.simulation as si
import pele_platform.Utilities.Parameters.parameters as pv


@dataclass
class GpcrLauncher:

    args: pv.ParametersBuilder

    def run_gpcr_simulation(self) -> pv.ParametersBuilder:
        # Set parameters for GPCR and launch simulation
        self._set_parameters()
        simulation_parameters = si.run_adaptive(self.args)
        return simulation_parameters

    def _set_parameters(self) -> None:
        # Set box and initial ligand position
        self.orthosteric_site = self.args.orthosteric_site
        self.initial_site = self.args.initial_site
        self.args.center_of_interface = self.initial_site
        self.args.randomize = True
