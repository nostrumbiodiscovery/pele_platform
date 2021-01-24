from pele_platform.Models.simulation_params_model import SimulationParamsModel
from pele_platform.Models.utils import PydanticProxy
from pele_platform.Models.yaml_parser_model import YamlParserModel
import pele_platform.Utilities.Helpers.helpers as hp
import pele_platform.constants.constants as cs

LOGFILE = '"simulationLogPath" : "$OUTPUT_PATH/logFile.txt",'

class SimulationParams(PydanticProxy):
    model_class = SimulationParamsModel

    def __init__(self, args: YamlParserModel):
        _args = args.model.dict()
        self.model_class.simulation_params = self.simulation_params
        self.initialize_model(_args)
        self.parse_args(args)

    def parse_args(self, args):
        self.optative_params(args)
        self.system_preparation_params(args)
        self.water_params(args)
        self.generate_inout_params(args)
        self.generate_pca_params(args)

    def generate_inout_params(self, args):
        # If User inputs exit condition or when doing an exit simulation
        if args.exit or args.in_out or args.in_out_soft:
            self.exit_value = (
                args.exit_value
                if args.exit_value
                else self.simulation_params.get("exit_value", 0.9)
            )
            self.exit_condition = (
                args.exit_condition
                if args.exit_condition
                else self.simulation_params.get("exit_condition", ">")
            )
            self.exit_trajnum = (
                args.exit_trajnum
                if args.exit_trajnum
                else self.simulation_params.get("exit_trajnum", 4)
            )
            self.bias_column = (
                args.bias_column
                if args.bias_column
                else self.simulation_params.get("bias_column", 6)
            )
            self.unbinding_block = cs.UNBINDING.format(
                self.bias_column,
                self.exit_value,
                self.exit_condition,
                self.exit_trajnum,
            )
            self.equilibration = "false"  # Lots of problems look into it
        else:
            self.unbinding_block = ""

    def generate_pca_params(self, args):
        if self.pca or self.pca_traj:
            self.anm_direction = (
                args.anm_direction
                if args.anm_direction
                else self.simulation_params.get("anm_direction", "oscillate")
            )
            self.anm_mix_modes = (
                args.anm_mix_modes
                if args.anm_mix_modes
                else self.simulation_params.get("anm_mix_modes", "doNotMixModes")
            )
            self.anm_picking_mode = (
                args.anm_picking_mode
                if args.anm_picking_mode
                else self.simulation_params.get("anm_picking_mode", "LOWEST_MODE")
            )
            self.anm_displacement = (
                args.anm_displacement
                if args.anm_displacement
                else self.simulation_params.get("anm_displacement", 2.0)
            )
            self.anm_modes_change = (
                args.anm_modes_change
                if args.anm_modes_change
                else self.simulation_params.get("anm_modes_change", 50)
            )
            self.anm_num_of_modes = (
                args.anm_num_of_modes
                if args.anm_num_of_modes
                else self.simulation_params.get("anm_num_of_modes", 1)
            )
            self.anm_relaxation_constr = (
                args.anm_relaxation_constr
                if args.anm_relaxation_constr
                else self.simulation_params.get("anm_relaxation_constr", 2)
            )
            self.remove_constraints = (
                args.remove_constraints
                if args.remove_constraints is not None
                else self.simulation_params.get("remove_constraints", True)
            )
            self.anm_freq = (
                args.anm_freq
                if args.anm_freq
                else self.simulation_params.get("anm_freq", 1)
            )

    def optative_params(self, args):
        self.equil_steps = (
            int(args.eq_steps / self.cpus) + 1
            if args.eq_steps
            else self.simulation_params.get("equilibration_steps", 1)
        )
        self.poses = (
            args.poses
            if args.poses
            else self.simulation_params.get("poses", self.cpus - 1)
        )

    def system_preparation_params(self, args):
        self.external_constraints = (
            hp.retrieve_constraints_for_pele(args.external_constraints, self.system)
            if args.external_constraints
            else []
        )

    def water_params(self, args):
        self.water_arg = (
            hp.retrieve_all_waters(self.system)
            if args.waters == "all_waters"
            else args.waters
        )  # IDs of waters
