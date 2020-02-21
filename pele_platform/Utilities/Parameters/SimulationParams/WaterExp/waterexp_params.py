import pele_platform.constants.pele_params as pcs


class WaterExp(object):

    def __init__(self, args):
        if args.water_exp:
            self.solvent = "VDGBNP"
            self.perturbation = ""
            self.perturbation_params(args)
            self.bias_column = args.bias_column if args.bias_column else self.simulation_params.get("bias_column", 4)
        if args.water_lig:
            if args.water_expl:
                self.parameters = pcs.WATER_LIG_EXPL
            else:
                self.parameters = pcs.WATER_LIG

