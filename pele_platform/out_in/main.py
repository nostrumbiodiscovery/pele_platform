from dataclasses import dataclass
import pele_platform.Adaptive.simulation as si
import pele_platform.Utilities.Helpers.helpers as hp
import pele_platform.Utilities.Parameters.parameters as pv
import pele_platform.Errors.custom_errors as ce


@dataclass
class OutInLauncher():

    args: pv.ParametersBuilder

    def run_gpcr_simulation(self) -> pv.ParametersBuilder:
        #Set parameters for gpcr and launch simulation
        self._check_mandatory_fields()
        self._set_parameters()
        simulation_parameters = si.run_adaptive(self.args)
        return simulation_parameters


    def _check_mandatory_fields(self):
        COMPULSORY_FLAGS = ["final_site", "initial_site"]
        for flag in COMPULSORY_FLAGS:
            if getattr(self.args, flag) is None:
                raise ce.OutInError(f"flag {flag} must be specified for out_in package")

        

    def _set_parameters(self) -> None:
        # Set box and initial ligand position
        self.final_site = self.args.final_site
        self.initial_site = self.args.initial_site
        self.args.center_of_interface = self.initial_site
        box_center, box_radius = hp.retrieve_box(
    self.args.system, self.initial_site, self.final_site,
    weights=[0.35, 0.65])
        self.args.box_center = self.args.box_center if self.args.box_center else box_center
        self.args.box_radius = self.args.box_radius if self.args.box_radius else box_radius
        self.args.randomize = True
