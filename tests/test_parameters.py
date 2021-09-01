import os

from pele_platform.Utilities.Parameters.parameters import Parameters
from pele_platform.Utilities.Helpers.yaml_parser import YamlParser
from pele_platform.constants import constants


def test_yaml_parser():

    yaml_file = os.path.join(constants.DIR, "Examples", "global", "input.yaml")
    breakpoint()
    parser_from_file = YamlParser.from_yaml(yaml_file)
    parser_from_file.read()

    params_dict = parser_from_file.to_dict()
    parser_from_dict = Parameters.from_dict(params_dict)
    parser_from_dict.read()

    assert parser_from_dict['global']
