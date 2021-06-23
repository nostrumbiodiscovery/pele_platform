import pele_platform.constants.constants as cs

class PCAParams(object):


    def generate_pca_params(self, args):
        self.pca = cs.PCA.format(args.pca) if args.pca else ""
        self.pca_traj = args.pca_traj
        if self.pca or self.pca_traj:
            self.anm_direction = args.anm_direction if args.anm_direction else self.simulation_params.get("anm_direction", "oscillate")
            self.anm_mix_modes = args.anm_mix_modes if args.anm_mix_modes else self.simulation_params.get("anm_mix_modes", "doNotMixModes")
            self.anm_picking_mode = args.anm_picking_mode if args.anm_picking_mode else self.simulation_params.get("anm_picking_mode", "LOWEST_MODE")
            self.anm_displacement = args.anm_displacement if args.anm_displacement else self.simulation_params.get("anm_displacement", 2.0)
            self.anm_modes_change = args.anm_modes_change if args.anm_modes_change else self.simulation_params.get("anm_modes_change", 50)
            self.anm_num_of_modes = args.anm_num_of_modes if args.anm_num_of_modes else self.simulation_params.get("anm_num_of_modes", 1)
            self.anm_relaxation_constr = args.anm_relaxation_constr if args.anm_relaxation_constr else self.simulation_params.get("anm_relaxation_constr", 2)
            self.remove_constraints = args.remove_constraints if args.remove_constraints is not None else self.simulation_params.get("remove_constraints", True)
            self.anm_freq = args.anm_freq if args.anm_freq else self.simulation_params.get("anm_freq", 1)
