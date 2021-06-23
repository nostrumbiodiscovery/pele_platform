import os
from dataclasses import dataclass
import subprocess
import pele_platform.Utilities.Helpers.helpers as hp
import pele_platform.Utilities.Parameters.parameters as pv
import pele_platform.Utilities.Helpers.plop_launcher as plop
import pele_platform.Errors.custom_errors as ce


@dataclass
class LigandParametrization:

    """
    Base class to generate ligand parameters
    1) Generate forcefield and rotamers from ligand
    2) Copy user's external files
    """

    env: pv.ParametersBuilder

    def generate(self) -> None:
        # Generate ligand forcefield parameters
        self.env.logger.info(
            "Creating template for residue {}".format(self.env.residue)
        )
        self.create_ligand_parameters()
        self.copy_external_parameters()
        self.env.logger.info(
            "Template {}z created\n\n".format(self.env.residue.lower())
        )

    def create_ligand_parameters(self) -> None:
        # create ligz template files
        with hp.cd(self.env.pele_dir):
            plop.parametrize_miss_residues(self.env)

    def copy_external_parameters(self) -> None:
        # copy user's files
        self.copy_ligand_forcefield_file()
        self.copy_ligand_rotamer_file()

    def copy_ligand_forcefield_file(self) -> None:
        # copy user's ligz template files
        for template_file in self.env.external_template:
            if not os.path.exists(template_file):
                raise ce.TemplateFileNotFound(f"File {template_file} not found")
            cmd_to_move_template = "cp {} {}".format(
                template_file, self.env.template_folder
            )
            subprocess.call(cmd_to_move_template.split())

    def copy_ligand_rotamer_file(self):
        # copy user's .rot.assign files
        for rotamer_file in self.env.external_rotamers:
            if not os.path.exists(rotamer_file):
                raise ce.RotamersFileNotFound(f"File {rotamer_file} not found")
            cmd_to_move_rotamer_file = "cp {} {}".format(
                rotamer_file, self.env.rotamers_folder
            )
            subprocess.call(cmd_to_move_rotamer_file.split())
