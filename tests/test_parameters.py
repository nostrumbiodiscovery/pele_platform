import os

from pele_platform.Utilities.Helpers.yaml_parser import YamlParser
from pele_platform.constants import constants


def test_yaml_parser():
    """
    Tests YamlParser methods:
     - to initialize the class from dict and YAML file
     - dump parameters as JSON and dict.
    """
    yaml_file = os.path.join(constants.DIR, "Examples", "global", "input.yaml")

    # Read from YAML
    parser_from_file = YamlParser.from_yaml(yaml_file)
    parser_from_file.read()

    # Dump to dict and JSON
    dictionary = parser_from_file.to_dict()
    assert dictionary["global"] is True  # assert it dumped as alias

    # Read from dictionary
    parser_from_dict = YamlParser.from_dict(dictionary)
    parser_from_dict.read()
    assert os.path.exists("input.yaml")  # ensure YAML got saved

    # Dump to JSON
    parser_from_dict.to_json()
