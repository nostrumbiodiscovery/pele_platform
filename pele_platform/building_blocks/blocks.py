from abc import ABC, abstractmethod
import os

from pele_platform.Utilities.Parameters import parameters


class Block(ABC):

    def __init__(
            self,
            parameters_builder: parameters.ParametersBuilder,
            options: dict,
            folder_name: str,
            env: parameters.Parameters = None,
    ) -> None:
        """
        Initialize Block class.

        Parameters
        ----------
        parameters_builder : ParametersBuilder
        options : dict
        folder_name : str
        env : Parameters
        """
        self.options = options
        self.folder_name = folder_name
        self.builder = parameters_builder
        self.original_dir = os.path.abspath(os.getcwd())
        self.env = env if env else self.builder.build_adaptive_variables(self.builder.initial_args)

    @abstractmethod
    def run(self):
        pass

    def set_working_folder(self):
        """
        Check if user specified a custom folder name for this simulation block. If not, use the automatically
        generated one.
        """
        if self.options:
            user_folder = self.options.get("working_folder", None)
            self.env.folder_name = user_folder if user_folder else self.folder_name
        else:
            self.env.folder_name = self.folder_name

        if self.env.pele_dir == self.builder.pele_dir:
            self.env.pele_dir = os.path.join(self.builder.pele_dir, self.env.folder_name)
        else:
            self.env.pele_dir = os.path.join(os.path.dirname(self.env.pele_dir), self.env.folder_name)

        self.env.inputs_dir = os.path.join(self.env.pele_dir, "input")
