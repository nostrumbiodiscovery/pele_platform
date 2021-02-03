from typing import Type, List
from pydantic import BaseModel, Field as PydanticField


def Field(
    *args,
    default=None,
    default_factory=None,
    value_from=None,
    value_from_simulation_params=None,
    simulation_params_default=None,
    can_be_falsy=False,
    tests_value=None,
    categories: List[str] = ["Other"],
    candidate_for_deprecation=False,
    **kwargs,
):
    """Extend signature for autocomplete."""
    kwargs = dict(
        **(
            {"default_factory": default_factory}
            if default_factory
            else {"default": default}
        ),
        value_from=value_from,
        value_from_simulation_params=value_from_simulation_params,
        simulation_params_default=simulation_params_default,
        can_be_falsy=can_be_falsy,
        tests_value=tests_value,
        categories=categories,
        candidate_for_deprecation=candidate_for_deprecation,
        **kwargs,
    )
    return PydanticField(*args, **kwargs)


class PydanticProxy:
    model_class: Type[BaseModel] = None
    model: BaseModel = None

    def initialize_model(self, data):
        self.model = self.model_class(**data)

    def __getattr__(self, key):
        if key in self.model.__fields__:
            return getattr(self.model, key)
        raise AttributeError

    def __setattr__(self, key, value):
        if self.model and key in self.model.__fields__:
            setattr(self.model, key, value)
        else:
            super().__setattr__(key, value)

    def __dir__(self):
        return dir(self.model)
