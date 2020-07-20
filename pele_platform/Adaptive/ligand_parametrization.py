from dataclasses import dataclass, field
from typing import List
import logging
import subprocess
import shutil
import offPELE.main as off
import pele_platform.Utilities.Helpers.helpers as hp
import pele_platform.Utilities.Parameters.pele_env as pv
import pele_platform.Utilities.Helpers.plop_launcher as plop
import pele_platform.constants.constants as cs

@dataclass
class LigandParametrization:

    '''
    Base class to generate ligand parameters
    1) Generate forcefield and rotamers from ligand
    2) Copy user's external files
    '''

    lig: str
    residue: str
    mae_lig: str=""
    external_templates: List = field(default_factory=list)
    external_rotamers: List = field(default_factory=list)
    pele_dir: str="."
    template_folder: str="."
    rotamers_folder: str="."
    core: int=-1
    mtor: int=4
    n: int=100000
    forcefield: str="OPLS2005"
    gridres: int=10
    logger: logging.Logger=None

    def generate(self) -> None:
        #Generate ligand forcefield parameters
        self.logger.info("Creating template for residue {}".format(self.residue))
        self.create_ligand_parameters()
        self.copy_external_parameters()
        self.logger.info("Template {}z created\n\n".format(self.residue.lower()))

    def create_ligand_parameters(self) -> None:
        #create ligz template files
        with hp.cd(self.pele_dir):
            if self.forcefield == cs.OPENFORCEFIELD:
                off.main(self.lig, as_datalocal=True)
            else:
                plop.parametrize_miss_residues(self)

    def copy_external_parameters(self) -> None:
        #copy user's files
        self.copy_ligand_forcefield_file()
        self.copy_ligand_rotamer_file()

    def copy_ligand_forcefield_file(self) -> None:
        #copy user's ligz template files
        for template_file in self.external_templates:
            cmd_to_move_template = "cp {} {}".format(template_file,  self.template_folder)
            subprocess.call(cmd_to_move_template.split())

    def copy_ligand_rotamer_file(self):
        #copy user's .rot.assign files
        for rotamer_file in self.external_rotamers:
            cmd_to_move_rotamer_file = "cp {} {}".format(rotamer_file, self.rotamers_folder)
            subprocess.call(cmd_to_move_rotamer_file.split())
