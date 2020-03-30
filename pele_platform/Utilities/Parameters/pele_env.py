import os
import shutil
import warnings
import random
import logging
import pele_platform.constants.constants as cs
import pele_platform.features.adaptive as adfs
import pele_platform.features.frag as frfs
import pele_platform.Utilities.Helpers.helpers as hp
from pele_platform.Utilities.Parameters.SimulationParams import simulation_params
from pele_platform.Utilities.Parameters.SimulationFolders import simulation_folders


class EnviroBuilder(simulation_params.SimulationParams, simulation_folders.SimulationPaths):
    """
        Base class wher the needed pele environment
        is build by creating folders and files
    """
    def __init__(self):
        pass

    def build_adaptive_variables(self, args):
        #DEFINE MAIN PATH
        pele_dir = os.path.abspath("{}_Pele".format(args.residue))
        if not args.folder:
            self.pele_dir = hp.is_repited(pele_dir) if args.restart in cs.FIRST_RESTART else hp.is_last(pele_dir)
            self.pele_dir = hp.is_repited(pele_dir) if not args.adaptive_restart else hp.is_last(pele_dir)
        else:
            self.pele_dir = os.path.abspath(args.folder)
        #####Define default variables, files and folder "HIDING VARIABLES " --> CHANGE#########
        for key, value in adfs.retrieve_software_settings(args, self.pele_dir).items():
            setattr(self, key, value)
        #####Initialize all variables by combining default and user input######
        simulation_params.SimulationParams.__init__(self, args)
        simulation_folders.SimulationPaths.__init__(self, args)
        for key, value in adfs.retrieve_software_settings(args, self.pele_dir).items():
            setattr(self, key, value)

    def build_frag_variables(self, args):
        #####Define default variables, files and folder "HIDING VARIABLES " --> CHANGE#########
        for key, value in frfs.retrieve_software_settings(args).items():
            setattr(self, key, value)
        #####Initialize all variables by combining default and user input######
        simulation_params.SimulationParams.__init__(self, args)
        for key, value in frfs.retrieve_software_settings(args).items():
            setattr(self, key, value)
        

    def create_files_and_folders(self):
        if not self.adaptive_restart:
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
