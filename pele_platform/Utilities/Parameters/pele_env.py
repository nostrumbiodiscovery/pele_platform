import os
import shutil
import warnings
import random
import logging
import pele_platform.constants as cs
import pele_platform.features as fs
import pele_platform.Utilities.Helpers.helpers as hp
from pele_platform.Utilities.Parameters.SimulationParams import simulation_params
from pele_platform.Utilities.Parameters.SimulationFolders import simulation_folders


class EnviroBuilder(simulation_params.SimulationParams, simulation_folders.SimulationPaths):
    """
        Base class wher the needed pele environment
        is build by creating folders and files
    """

    def __init__(self, args):
        self.build_variables(args)

    def build_variables(self, args):
        simulation_params.SimulationParams.__init__(self, args)
        simulation_folders.SimulationPaths.__init__(self, args)
        #####Define files and folders "HIDING VARIABLES " --> CHANGE#########
        for key, value in fs.retrieve_software_settings(args, self.pele_dir).items():
            setattr(self, key, value)

    @classmethod
    def build_env(cls, args):
        env = cls(args)
        env.create()
        return env

    def create(self):
        if self.restart in cs.FIRST_RESTART:
            self.create_folders()
            self.create_files()
            self.create_logger()
        else:
            self.create_logger()



    def create_folders(self):
        """
            Create pele folders
        """

        for folder in self.folders:
            hp.create_dir(self.pele_dir, folder)

    def create_files(self):
        """
            Copy templates
        """

        # Actions
        for file, destination_name in zip(self.files, self.file_names):
            shutil.copy(file, os.path.join(self.pele_dir, destination_name))


    def create_logger(self):
        log_name = os.path.join(self.pele_dir, "{}.log".format(self.residue))
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)
        formatter = logging.Formatter("%(asctime)s:%(levelname)s:%(message)s")
        if self.restart in ["all", "glide" ]:
            file_handler = logging.FileHandler(log_name, mode='w')
        else:
            file_handler = logging.FileHandler(log_name, mode='a')
        file_handler.setFormatter(formatter)
        self.logger.addHandler(file_handler)
