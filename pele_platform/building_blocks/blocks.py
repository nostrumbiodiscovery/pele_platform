from abc import ABC, abstractmethod
import os

from pele_platform.context import context


class Block(ABC):

    def __init__(
            self,
            options: dict,
            folder_name: str,
    ) -> None:
        """
        Initialize Block class.

        Parameters
        ----------
        options : dict
            Dictionary with arguments to override default Block parameters.
        folder_name : str
            Name of the folder in which the Block will be executed (one level under LIG_Pele).
        """
        self.options = options
        self.folder_name = folder_name
        self.original_dir = os.path.abspath(os.getcwd())
        if not context.parameters:
            context.parameters_builder.build_adaptive_variables()

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
            context.parameters.folder_name = user_folder if user_folder else self.folder_name
        else:
            context.parameters.folder_name = self.folder_name

        if context.parameters.pele_dir == context.parameters_builder.pele_dir:
            context.parameters.pele_dir = os.path.join(context.parameters_builder.pele_dir,
                                                       context.parameters.folder_name)
        else:
            context.parameters.pele_dir = os.path.join(os.path.dirname(context.parameters.pele_dir),
                                                       context.parameters.folder_name)

        context.parameters.inputs_dir = os.path.join(context.parameters.pele_dir, "input")
