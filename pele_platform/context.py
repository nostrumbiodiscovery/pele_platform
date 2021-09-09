import copy
from dataclasses import dataclass

if __name__ == "__main__":
    from pele_platform.Utilities.Helpers.yaml_parser import YamlParser
    from pele_platform.Utilities.Parameters.parameters import Parameters, ParametersBuilder


@dataclass
class Context:

    yaml_parser: "YamlParser" = None
    parameters: "Parameters" = None
    parameters_builder: "ParametersBuilder" = None

    def reset(self):
        """
        Resets Context instance.
        """
        self.parameters = None
        self.yaml_parser = None
        self.parameters_builder = None

    @property
    def parameters_copy(self):
        return copy.deepcopy(self.parameters)


context = Context()
