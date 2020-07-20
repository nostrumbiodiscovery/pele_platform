import pele_platform.Adaptive.solvent as sv


class FragSolvent():

    def __init__(self, args):
        self.solvent = "VDGBNP"
        self.solvent = sv.ImplicitSolvent(self.solvent, self.obc_tmp,
             self.template_folder, self.obc_file, self.logger,
             self.forcefield)
        self.solvent.generate()

