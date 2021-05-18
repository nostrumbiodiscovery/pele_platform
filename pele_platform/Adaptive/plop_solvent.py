import os
from dataclasses import dataclass
import glob
import shutil
import logging
import pele_platform.Utilities.Helpers.solventOBCParamsGenerator as obc


@dataclass
class ImplicitSolvent:

    solvent: str
    obc_tmp: str
    template_folder: str
    obc_file: str
    logger: logging.Logger

    def generate(self):
        self.logger.info("Setting implicit solvent: {}".format(self.solvent))
        self.set_implicit_solvent()
        self.logger.info("Implicit solvent set\n\n".format(self.solvent))

    def set_implicit_solvent(self):
        if self.solvent == "OBC":  # OBC from ...
            shutil.copy(self.obc_tmp, self.obc_file)
            for template in glob.glob(os.path.join(self.template_folder, "*")):
                obc.main(template, self.obc_file)
        else:
            pass  # SVG by default
