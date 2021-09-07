from difflib import SequenceMatcher
import yaml

from pele_platform.features.adaptive import SOFTWARE_CONSTANTS
from pele_platform.Errors import custom_errors
from pele_platform.Models.utils import PydanticProxy
from pele_platform.Models.yaml_parser_model import YamlParserModel


def _yaml_error_wrapper(error):
    """
    Wraps YAML errors into a more human-friendly format, making customs suggestions about potential issues.
    This should be expanded in the future when more issues get reported by the users.
    """
    custom_yaml_errors = {
        "expected '<document start>', but found '<block mapping start>'": "Please ensure every key in input.yaml is "
        "followed by a colon and "
        "a space. There seem to be some issues on line {}, character {}.",
        "found character '\\t' that cannot start any token": "Please remove any trailing tabs from input.yaml, there "
        "seem to be one on line {}, character {}.",
    }

    custom = custom_yaml_errors.get(str(error.problem).strip(), None)

    if custom:
        line = error.problem_mark.line + 1
        character = int(error.problem_mark.column) + 1
        raise yaml.YAMLError(custom.format(line, character))
    else:
        raise error


class YamlParser(PydanticProxy):

    def __init__(self, yaml_file=None, dictionary=None):
        self.yaml_file = yaml_file
        self.dictionary = dictionary
        self.model_class = YamlParserModel
        self.data = None
        self.valid_flags = None

    def read(self):
        """
        Parses and checks data provided by the user (either in the form of YAML file or dictionary).
        """
        self.data = self.dictionary if self.dictionary else self._parse_yaml()
        self.initialize_model(self.data)
        self.valid_flags = set(
            field.alias if field.alias != field.name else field.name
            for field in self.model.__fields__.values()
        )
        self._check(self.data)
        self._check_multiple_simulations()

    @classmethod
    def from_yaml(cls, user_file):
        """
        Initializes YamlParser class from YAML file provided by the user.

        Parameters
        ----------
        user_file : str
            Path to YAML file.

        Returns
        -------
            YamlParser class.
        """
        return YamlParser(yaml_file=user_file)

    @classmethod
    def from_dict(cls, user_dict):
        """
        Initializes YamlParser class from a dictionary of parameters and saves parameters to input.yaml file.

        Parameters
        ----------
        user_dict : dict
            Dictionary with simulation parameters.

        Returns
        -------
            YamlParser class.
        """
        parser = YamlParser(dictionary=user_dict)
        parser.yaml_file = "input.yaml"

        with open(parser.yaml_file, "w+") as file:
            yaml.dump(user_dict, file)

        return parser

    def to_dict(self):
        """
        Dumps all YamlParser parameters to a dictionary.

        Returns
        -------
            YamlParser parameters as dict.
        """
        return self.model.dict(exclude_none=True, by_alias=True)

    def to_json(self):
        """
        Dumps all YamlParser parameters to JSON.

        Returns
        -------
            YamlParser parameters in JSON format.
        """
        return self.model.json(exclude_none=True, by_alias=True)

    def initialize_model(self, data):
        """
        Initializes pydantic model.

        Parameters
        ----------
        data : dict
            Dictionary of user-defined parameters.
        """
        if not isinstance(data, dict):
            raise custom_errors.WrongYamlFile(f"Input file: {self.yaml_file} does not look like a correct YAML file.")

        super().initialize_model(data)

    def _parse_yaml(self) -> dict:
        """
        Converts YAML file to dictionary.

        Returns
        -------
            Dictionary of user-defined parameters.
        """
        with open(self.yaml_file, "r") as stream:
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
            raise custom_errors.MultipleSimulationTypes(
                f"You cannot select multiple simulation types in input.yaml, please select one of "
                f"{', '.join(specified_simulations)}."
            )

    def _check(self, data) -> None:
        """
        Checks if all user-defined flags are valid.
        Parameters
        ----------
        data : dict
            Dictionary of user-defined data.
        Raises
        ------
        KeyError when invalid flag is identified.
        """
        for key in data.keys():
            if key not in self.valid_flags:
                raise KeyError(self._recommend(key))

    def _recommend(self, key):
        """
        Recommends the most similar flag to replace the invalid one.

        Parameters
        ----------
        key : str
            Invalid flag.

        Returns
        -------
            String with a suggestion to fix the invalid flag.
        """
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


class MostSimilarFlag:

    def __init__(self, name):
        self.name = name
        self.distance = None

    def calculate_distance(self, key):
        self.distance = SequenceMatcher(None, self.name, key).ratio()
