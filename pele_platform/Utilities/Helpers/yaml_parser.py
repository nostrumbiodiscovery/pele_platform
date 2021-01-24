from dataclasses import dataclass
from difflib import SequenceMatcher
from pele_platform.Errors import custom_errors
from pele_platform.Models.utils import PydanticProxy
from pele_platform.Models.yaml_parser_model import YamlParserModel
import yaml

@dataclass
class YamlParser(PydanticProxy):
    yamlfile: str
    valid_flags: set = None
    model_class = YamlParserModel

    def read(self) -> None:
        data = self._parse_yaml()
        self.initialize_model(data)
        self.valid_flags = set(
            field.alias if field.alias != field.name else field.name
            for field in self.model.__fields__.values()
        )
        self._check(data)

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
            except yaml.YAMLError as exc:
                raise (exc)
        return data

    def _check(self, data) -> None:
        # Check if valids in yaml file are valids
        for key in data.keys():
            if key not in self.valid_flags:
                raise KeyError(self._recommend(key))

    def _recommend(self, key):
        most_similar_flag = None
        for valid_key in self.valid_flags:
            flag = Most_Similar_Flag(valid_key)
            flag.calculate_distance(key)
            if not most_similar_flag:
                most_similar_flag = flag
            else:
                if flag.distance > most_similar_flag.distance:
                    most_similar_flag = flag
        return f"Incorrect flag {key}. Did you mean {most_similar_flag.name}?"


@dataclass
class Most_Similar_Flag:

    name: str

    def calculate_distance(self, key):
        self.distance = SequenceMatcher(None, self.name, key).ratio()
