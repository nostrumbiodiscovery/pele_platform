from dataclasses import dataclass
import pele_platform.Adaptive.simulation as si
import pele_platform.Utilities.Helpers.helpers as hp
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
        box_center, box_radius = hp.retrieve_box(
    self.args.system, self.initial_site, self.orthosteric_site,
    weights=[0.35, 0.65])
        self.args.box_center = self.args.box_center if self.args.box_center else box_center
        self.args.box_radius = self.args.box_radius if self.args.box_radius else box_radius
        self.args.randomize = True
