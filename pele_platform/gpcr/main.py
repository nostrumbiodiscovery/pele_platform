from dataclasses import dataclass
import pele_platform.Adaptive.simulation as si
import pele_platform.Utilities.Helpers.helpers as hp
import pele_platform.Utilities.Parameters.pele_env as pv


@dataclass
class GpcrLauncher():

    args: pv.EnviroBuilder 

    def run_gpcr_simulation(self) -> pv.EnviroBuilder:
        #Set parameters for gpcr and launch simulation
        self._set_parameters()
        simulation_parameters = si.run_adaptive(self.args)
        return simulation_parameters
        

    def _set_parameters(self) -> None:
        # Set box and initial ligand position
        self.orthosteric_site = self.args.orthosteric_site
        self.initial_site = self.args.initial_site
        self.args.center_of_interface = self.initial_site
        self.args.box_center, self.args.box_radius = hp.retrieve_box(
    self.args.system, self.initial_site, self.orthosteric_site,
    weights=[0.35, 0.65])
        self.args.randomize = True
