import os
import shutil
import logging
from pele_platform.constants import constants
from pele_platform.features import adaptive
from pele_platform.features import frag
from pele_platform.Utilities.Helpers import helpers
from pele_platform.Utilities.Parameters.SimulationParams import simulation_params
from pele_platform.Utilities.Parameters.SimulationFolders import simulation_folders


class ParametersBuilder(object):
    """
    It builds a Parameters instance and it creates the required folders and
    files

    .. todo::
       * Consider removing the functionality about the creation of files
         and folders outside this method.
    """
    def __init__(self):
        """
        It initializes a ParametersBuilder object.
        """
        self.parameters = None
        self._initialized = False

    def build_adaptive_variables(self, args):
        """
        It builds the parameters for adaptive, according to the arguments
        that are supplied, and returns the corresponding Parameters
        instance.

        Parameters
        ----------
        args :
        """
        # Define main path
        pele_dir = os.path.abspath("{}_Pele".format(args.residue))
        if not args.folder:
            self.pele_dir = helpers.is_repeated(pele_dir) if args.restart in constants.FIRST_RESTART else helpers.is_last(pele_dir)
            self.pele_dir = helpers.is_repeated(pele_dir) if not args.adaptive_restart else helpers.is_last(pele_dir)
        else:
            self.pele_dir = os.path.abspath(args.folder)
        #####Define default variables, files and folder "HIDING VARIABLES " --> CHANGE#########
        for key, value in adaptive.retrieve_software_settings(args, self.pele_dir).items():
            setattr(self, key, value)
        #####Initialize all variables by combining default and user input######
        simulation_params.SimulationParams.__init__(self, args)
        simulation_folders.SimulationPaths.__init__(self, args)
        for key, value in adaptive.retrieve_software_settings(args, self.pele_dir).items():
            setattr(self, key, value)

    def build_frag_variables(self, args):
        # Define default variables, files and folder "HIDING VARIABLES " --> CHANGE
        for key, value in frag.retrieve_software_settings(args).items():
            setattr(self, key, value)
        # Initialize all variables by combining default and user input
        simulation_params.SimulationParams.__init__(self, args)
        for key, value in frag.retrieve_software_settings(args).items():
            setattr(self, key, value)
        self.create_logger(".")

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
            helpers.create_dir(self.pele_dir, folder)

    def create_files(self):
        """
            Copy templates
        """

        # Actions
        for file, destination_name in zip(self.files, self.file_names):
            shutil.copy(file, os.path.join(self.pele_dir, destination_name))

    def create_logger(self, directory=None):
        directory = directory if directory else self.pele_dir
        log_name = os.path.join(directory, "{}.log".format(self.residue))
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)
        self.logger.propagate = False
        formatter = logging.Formatter("%(asctime)s:%(levelname)s:%(message)s")
        if not self.adaptive_restart:
            file_handler = logging.FileHandler(log_name, mode='w')
        else:
            file_handler = logging.FileHandler(log_name, mode='a')
        file_handler.setFormatter(formatter)
        self.logger.addHandler(file_handler)

    @property
    def initialized(self):
        """
        It returns the initialization state of the builder. If this builder
        has already built a Parameters object before, it will return a True.
        Otherwise, it will return a False.

        Returns
        -------
        initialized : bool
            Whether this ParametersBuilder has already been used to build
            a Parameters object or not
        """
        return self._initialized


class Parameters(simulation_params.SimulationParams,
                 simulation_folders.SimulationPaths):
    """
    Base parameters class where the PELE Platform parameters are stored
    and manipulated.

    .. todo::
       * This class should be serializable. For example, we should be able to
         represent it as a json format at any time. We need to implement
         an abstract class for this purpose, something like in:
         https://github.com/openforcefield/openff-toolkit/blob/master/openff/toolkit/utils/serialization.py#L25
    """

    def __init__(self, args):
        """
        It initializes a Parameters object.

        Parameters
        ----------
        args : The
        """
        # Initialize the parameters from parent classes
        super().__init__()

        # CA interval is initialized to None
        self.ca_interval = None