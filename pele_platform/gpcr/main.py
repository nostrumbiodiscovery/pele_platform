import pele_platform.Adaptive.simulation as si
import pele_platform.Utilities.Helpers.helpers as hp



class GpcrLauncher():

    def __init__(self, args):
        self.orthosteric_site = args.orthosteric_site
        self.initial_site = args.initial_site
        self.args = args

    def run_gpcr_simulation(self):
        self._set_parameters()
        simulation_parameters = si.run_adaptive(self.args)
        return simulation_parameters
        

    def _set_parameters(self):
        #Set box and system initial parameters
        self.args.center_of_interface = self.initial_site
        self.args.box_center, self.args.box_radius = hp.retrieve_box(
    self.args.system, self.initial_site, self.orthosteric_site,
    weights=[0.35, 0.65])
        self.args.randomize = True
        print(self.args.box_center, self.args.box_radius)
