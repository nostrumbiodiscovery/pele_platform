import copy
from dataclasses import dataclass

if __name__ == "__main__":  # dirty... find a better solution
    from pele_platform.Utilities.Helpers.yaml_parser import YamlParser
    from pele_platform.Utilities.Parameters.parameters import Parameters, ParametersBuilder


@dataclass
class Context:

    yaml_parser: "YamlParser" = None
    parameters: "Parameters" = None
    parameters_builder: "ParametersBuilder" = None

    def reset_parameters(self):
        """
        Resets Context instance.
        """
        self.parameters = None
        self.yaml_parser = None
        self.parameters_builder = None

        # TODO: Context manager with new_context...

    @property
    def parameters_copy(self):
        return copy.deepcopy(self.parameters)

    def json(self):
        """
        Dumps Parameters to JSON.
        """
        pass

    def dict(self):
        """
        Dumps parameters to a dictionary.
        """
        pass

    def build_parameters(self, yaml_file=None, user_dict=None):
        # TODO: get rid of those imports
        from pele_platform.Utilities.Helpers.yaml_parser import YamlParser
        from pele_platform.Utilities.Parameters.parameters import ParametersBuilder

        if yaml_file:
            self.yaml_parser = YamlParser.from_yaml(yaml_file)
        else:
            self.yaml_parser = YamlParser.from_dict(user_dict)

        self.yaml_parser.read()
        self.parameters_builder = ParametersBuilder()


context = Context()
