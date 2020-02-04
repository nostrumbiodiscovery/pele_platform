import pele_platform.constants.constants as cs

class PCAParams(object):


    def __init__(self, args):
        self.pca = cs.PCA.format(args.pca) if args.pca else ""
        self.pca_traj = args.pca_traj
        if self.pca or self.pca_traj:
            self.anm_direction = args.anm_direction if args.anm_direction else "oscillate"
            self.anm_mix_modes = args.anm_mix_modes if args.anm_mix_modes else "doNotMixModes"
            self.anm_picking_mode = args.anm_picking_mode if args.anm_picking_mode else "LOWEST_MODE"
            self.anm_displacement = args.anm_displacement if args.anm_displacement else 1.75
            self.anm_modes_change = args.anm_modes_change if args.anm_modes_change else 8
            self.anm_num_of_modes = args.anm_num_of_modes if args.anm_num_of_modes else 1
            self.anm_relaxation_constr = args.anm_relaxation_constr if args.anm_relaxation_constr else 2
            self.remove_constraints = args.remove_constraints if not args.remove_constraints else True
        else:
            self.anm_direction = args.anm_direction if args.anm_direction else "random"
            self.anm_mix_modes = args.anm_mix_modes if args.anm_mix_modes else "mixMainModeWithOthersModes"
            self.anm_picking_mode = args.anm_picking_mode if args.anm_picking_mode else "RANDOM_MODE"
            self.anm_displacement = args.anm_displacement if args.anm_displacement else 0.75
            self.anm_modes_change = args.anm_modes_change if args.anm_modes_change else 4
            self.anm_num_of_modes = args.anm_num_of_modes if args.anm_num_of_modes else 6
            self.anm_relaxation_constr = args.anm_relaxation_constr if args.anm_relaxation_constr else 0.5
            self.remove_constraints = args.remove_constraints if args.remove_constraints else False
