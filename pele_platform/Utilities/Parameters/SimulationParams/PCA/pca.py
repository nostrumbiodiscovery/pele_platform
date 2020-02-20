import pele_platform.constants.constants as cs

class PCAParams(object):


    def __init__(self, args):
        self.pca = cs.PCA.format(args.pca) if args.pca else ""
        self.pca_traj = args.pca_traj
        if self.pca or self.pca_traj:
            self.anm_direction = self.simulation_params.get("anm_direction", "oscillate")
            self.anm_mix_modes = self.simulation_params.get("anm_mix_modes", "doNotMixModes")
            self.anm_picking_mode = self.simulation_params.get("anm_picking_mode", "LOWEST_MODE")
            self.anm_displacement = self.simulation_params.get("anm_displacement", 1.75)
            self.anm_modes_change = self.simulation_params.get("anm_modes_change", 8)
            self.anm_num_of_modes = self.simulation_params.get("anm_num_of_modes", 1)
            self.anm_relaxation_constr = self.simulation_params.get("anm_relaxation_constr", 2)
            self.remove_constraints = self.simulation_params.get("remove_constraints", True)
