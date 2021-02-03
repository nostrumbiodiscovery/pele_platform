from collections import defaultdict
from pele_platform.Utilities.Helpers.yaml_parser import YamlParserModel
from pele_platform.Utilities.Parameters.SimulationParams.simulation_params import (
    SimulationParamsModel,
)
try:
    from jinja2 import Template
except ImportError:

    pass
import os


def get_template(template):
    return Template(open(os.path.join(os.path.dirname(__file__), template)).read())


def title(value):
    value = value.replace("_", " ")
    if value == value.lower():
        return value.title()
    return value


def parse_field(field):
    extra = field.field_info.extra
    vfsp = extra.get("value_from_simulation_params")
    return {
        "categories": extra.get("categories"),
        "title": field.field_info.title
        if field.field_info.title
        else (
            title(field.alias)
            + (f" / {title(field.name)}" if field.name != field.alias else "")
        ),
        "name": field.name,
        "alias": field.alias,
        "description": field.field_info.description,
        "default": repr(field.default) if field.default is not None else f"{field.default_factory.__name__}()" if field.default_factory is not None else None,
        "type_": field.type_.__name__ if type(field.type_) == type else str(field.type_).replace('typing.',''),
        **extra,
        "tests_value": repr(extra.get("tests_value"))
        if extra.get("tests_value") is not None
        else None,
        "simulation_params_default": repr(extra.get("simulation_params_default"))
        if extra.get("simulation_params_default") is not None
        else None,
        "value_from_simulation_params": (field.name if vfsp is True else vfsp)
        if vfsp
        else None,
        "field": field,
        "validators": [
            f"``{v[0]}`` from `{v[1]}`"
            for v in [
                validator.__qualname__.split(".")[::-1] + ["Pydantic"]
                for validator in field.validators
                + (field.post_validators if field.post_validators else [])
            ]
        ],
        "from_simulation_params": field.name not in YamlParserModel.__fields__,
    }


def render_yaml_parser():
    fields = [parse_field(field) for field in SimulationParamsModel.__fields__.values()]

    categories = defaultdict(list)
    for field in fields:
        for category in field['categories']:
            categories[category].append(field)

    template = get_template("./docgen.rst")
    return template.render(categories=categories)
