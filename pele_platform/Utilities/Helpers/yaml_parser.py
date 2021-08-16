from dataclasses import dataclass
from pele_platform.Errors.custom_errors import (
    LigandNameNotSupported,
    MultipleSimulationTypes,
)
from pele_platform.features.adaptive import SOFTWARE_CONSTANTS
from difflib import SequenceMatcher
from pele_platform.Errors import custom_errors
from pele_platform.Models.utils import PydanticProxy
from pele_platform.Models.yaml_parser_model import YamlParserModel
import yaml


def _yaml_error_wrapper(error):
    """
    Wraps YAML errors into a more human-friendly format, making customs suggestions about potential issues.
    This should be expanded in the future when more issues get reported by the users.
    """
    custom_errors = {
        "expected '<document start>', but found '<block mapping start>'": "Please ensure every key in input.yaml is "
        "followed by a colon and "
        "a space. There seem to be some issues on line {}, character {}.",
        "found character '\\t' that cannot start any token": "Please remove any trailing tabs from input.yaml, there "
        "seem to be one on line {}, character {}.",
    }

    custom = custom_errors.get(str(error.problem).strip(), None)

    if custom:
        line = error.problem_mark.line + 1
        character = int(error.problem_mark.column) + 1
        raise yaml.YAMLError(custom.format(line, character))
    else:
        raise error


@dataclass
class YamlParser(PydanticProxy):
    yamlfile: str
    valid_flags: set = None
    model_class = YamlParserModel
    data: dict = None

    def read(self) -> None:
        self.data = self._parse_yaml()
        self.initialize_model(self.data)
        self.valid_flags = set(
            field.alias if field.alias != field.name else field.name
            for field in self.model.__fields__.values()
        )
        self._check(self.data)
        self._check_multiple_simulations()

    def initialize_model(self, data):
        if not isinstance(data, dict):
            raise custom_errors.WrongYamlFile(
                "Input file: {} does not look like a correct yml file".format(
                    self.yamlfile
                )
            )
        super().initialize_model(data)

    def _parse_yaml(self) -> dict:
        # Retrieve raw info from yaml
        with open(self.yamlfile, "r") as stream:
            try:
                data = yaml.safe_load(stream)
            except Exception as error:
                _yaml_error_wrapper(error)
        return data

    def _check_multiple_simulations(self):
        """
        Checks if the user specified more than one simulation type in YAML.

        Raises
        -------
        MultipleSimulationTypes if more than one simulation type set in YAML.
        """
        available_simulations = SOFTWARE_CONSTANTS.get("simulation_params", {})
        specified_simulations = [
            key for key in self.data.keys() if key in available_simulations.keys()
        ]

        if len(specified_simulations) > 1:
            raise MultipleSimulationTypes(
                f"You cannot select multiple simulation types in input.yaml, please select one of "
                f"{', '.join(specified_simulations)}."
            )

    def _check(self, data) -> None:
        # Check if valids in yaml file are valids
        for key in data.keys():
            if key not in self.valid_flags:
                raise KeyError(self._recommend(key))

    def _recommend(self, key):
        most_similar_flag = None
        for valid_key in self.valid_flags:
            flag = MostSimilarFlag(valid_key)
            flag.calculate_distance(key)
            if not most_similar_flag:
                most_similar_flag = flag
            else:
                if flag.distance > most_similar_flag.distance:
                    most_similar_flag = flag
        return f"Incorrect flag {key}. Did you mean {most_similar_flag.name}?"


@dataclass
class MostSimilarFlag:
    name: str

    def calculate_distance(self, key):
        self.distance = SequenceMatcher(None, self.name, key).ratio()
