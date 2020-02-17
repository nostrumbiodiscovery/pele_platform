import pele_platform.constants.pele_params as pcs


class WaterExp(object):

    def __init__(self, args):
        if args.water_exp:
            self.solvent = "VDGBNP"
            self.perturbation = ""
            self.perturbation_params(args)
        if args.water_lig:
            if args.water_expl:
                self.parameters = pcs.WATER_LIG_EXPL
            else:
                self.parameters = pcs.WATER_LIG

