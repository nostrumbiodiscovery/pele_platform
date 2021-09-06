from pele_platform.Models.simulation_params_model import SimulationParamsModel
from pele_platform.Models.utils import PydanticProxy
from pele_platform.Models.yaml_parser_model import YamlParserModel


class SimulationParams(PydanticProxy):
    model_class = SimulationParamsModel

    def __init__(self, args: YamlParserModel):
        _args = args.model.dict()
        self.model_class.simulation_params = self.simulation_params
        self.initialize_model(_args)
        self.parse_args(args)

    def to_dict(self):
        """
        Dumps all SimulationParams to a dictionary.

        Returns
        -------
            SimulationParams parameters as dict.
        """
        return self.model.dict()

    def to_json(self):
        """
        Dumps all SimulationParams parameters to JSON.

        Returns
        -------
            SimulationParams parameters in JSON format.
        """
        return self.model.json()
