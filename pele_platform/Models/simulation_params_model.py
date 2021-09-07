import os
import glob
from typing import Any, List
from pydantic import validator, root_validator

from pele_platform.constants import constants
import pele_platform.Utilities.Helpers.helpers as hp
from pele_platform.Models.utils import Field
from pele_platform.Models.yaml_parser_model import YamlParserModel
from pele_platform.features import adaptive as features


class SimulationParamsModel(YamlParserModel):
    frag_pele: bool = Field()
    complexes: str = Field()
    frag_pele_steps: int = Field()
    output_path: str = Field()
    water: str = Field()
    ligand: str = Field(default=constants.LIGAND)
    spython: str = Field(
        default=os.path.join(constants.SCHRODINGER, "utilities/python")
    )
    lig: str = Field(value_from="mae_lig")
    sasa_max: Any = Field()
    sasa_min: Any = Field()
    clusters: int = Field(alias="clust", tests_value=2)
    allow_empty_selectors: bool = Field(value_from="water_empty_selector")
    templates: Any = Field(
        default=os.path.abspath(
            os.path.join(os.path.dirname(os.path.dirname(__file__)), "PeleTemplates")
        )
    )
    simulation_type: str = Field(default="pele")
    xtc: bool = Field()
    pdb: bool = Field()
    constraints: List[str] = Field()
    water_energy: Any = Field()
    sidechain_perturbation: bool = Field(default=False)
    met_interaction_restrictions: str = Field()
    covalent_sasa: str = Field()
    max_trials_for_one: int = Field()
    conformation_perturbation: str = Field()
    equilibration_mode: str = Field()
    water_ids_to_track: List[str] = Field(default=list())
    inputs_dir: str = Field()
    residue_type: str = Field()
    inputs: List[str] = Field()
    ligand_ref: str = Field()
    mpi_params: str = Field()
    logfile: str = Field()
    sasa: str = Field(
        value_from_simulation_params=True,
        simulation_params_default=constants.SASA,
    )
    parameters: str = Field(
        value_from_simulation_params="params",
        simulation_params_default=True,
    )
    selection_to_perturb: str = Field(
        value_from_simulation_params=True,
        simulation_params_default=constants.SELECTION_TO_PERTURB,
    )
    unbinding_block: str = Field()
    water_arg: str = Field()

    @validator("*", pre=True, always=True)
    def set_value_from_simulation_parameters(cls, v, field):

        extra = field.field_info.extra
        can_be_falsy = extra.get("can_be_falsy")

        if (can_be_falsy and v is None) or (not can_be_falsy and not v):
            value_from_simulation_params = extra.get("value_from_simulation_params")

            if value_from_simulation_params is True:  # TODO: IS?
                value_from_simulation_params = field.name

            if value_from_simulation_params:
                value = cls.simulation_params.get(
                    value_from_simulation_params, extra.get("simulation_params_default")
                )
                if value == []:
                    return []  # Prevent mutation bugs
                return value
        return v

    @validator("adaptive", always=True)
    def set_adaptive(cls, v, values):
        if values.get("package") in features.all_simulations:
            return True
        return False

    @root_validator
    def set_frag_pele_parameters(cls, values):
        """
        Update all placeholders in the template in case we're running FragPELE instead of AdaptivePELE.
        """
        if values.get("frag_core"):
            values["frag_pele"] = True
            values["complexes"] = "$PDB"
            values["frag_pele_steps"] = "$STEPS"
            values["output_path"] = "$RESULTS_PATH"
        else:
            values["frag_pele"] = False
            values["complexes"] = "$COMPLEXES"
            values["frag_pele_steps"] = "$PELE_STEPS"
            values["output_path"] = "$OUTPUT_PATH"

        return values

    @root_validator
    def set_out_in_parameters(cls, values):

        if values["exit"] or values["in_out"] or values["in_out_soft"]:
            values["equilibration"] = False
            values["unbinding_block"] = constants.UNBINDING.format(
                values["bias_column"],
                values["exit_value"],
                values["exit_condition"],
                values["exit_trajnum"],
            )
        else:
            values["unbinding_block"] = ""

        return values

    @validator("input", always=True)
    def set_input_glob(cls, v, values):
        """
        Extract all input PDBs, if 'system' contains an asterisk.

        Parameters
        -----------
        v : str
            Value set for argument being checked.
        values : dict
            Dictionary of all parser parameters.

        Returns
        -------
            List of paths to PELE inputs.
        """
        system = values.get("system")
        if system and "*" in system:
            return glob.glob(system)
        return v

    @validator("system", always=True)
    def set_system_glob(cls, v):
        """
        Set a single PDB file as system, in case 'system' contains an asterisk.

        Parameters
        -----------
        v : str
            Value set for argument being checked.

        Returns
        -------
            Path to the first PELE input.
        """
        if v and "*" in v:
            return glob.glob(v)[0]
        return v

    @validator("system", "residue", "chain", always=True)
    def validate_adaptive_required_fields(cls, v, values, field):
        """
        Ensure that fields required for adaptive simulation were set.

        Parameters
        -----------
        v : str
            Value set for argument being checked.
        values : dict
            Dictionary of all parser parameters.
        field : Field
            Object describing parameters of the YAML argument.

        Raises
        -------
            AssertionError if one of the parameters is missing.

        Returns
        --------
            Value of the argument being checked, if it's correct.
        """
        if values.get("adaptive") and not v:
            raise AssertionError(f"User must define {field.name} to use AdaptivePELE.")
        return v

    @validator("perturbation", always=True)
    def set_perturbation(cls, v):
        if v is False:
            return ""
        return cls.simulation_params.get("perturbation", constants.PERTURBATION)

    @validator("sidechain_perturbation", always=True)
    def set_sidechain_perturbation(cls, v, values):
        if v or values.get("covalent_residue"):
            return True
        return False

    @validator("proximityDetection", always=True)
    def set_proximityDetection(cls, v):
        if v is False:
            return "false"
        return cls.simulation_params.get("proximityDetection", True)

    @validator(
        "selection_to_perturb",
        "parameters",
        "ligand",
        "binding_energy",
        "sasa",
        always=True,
    )
    def only_with_perturbation(cls, v, values):
        if not values.get("perturbation"):
            return ""
        return v

    @validator("conformation_freq", always=True)
    def only_with_conformation_perturbation(cls, v, values):
        if values.get("ligand_conformations"):
            return constants.CONFORMATION_FREQUENCY.format(v)
        return ""

    @validator("spython", always=True)
    def check_spython_path(cls, v):
        if not os.path.exists(v):
            return os.path.join(constants.SCHRODINGER, "run")
        return v

    @validator("lig", always=True)
    def set_lig(cls, v, values):
        if v:
            return v
        return "{}.mae".format(values.get("residue"))

    @validator("iterations", always=True)
    def set_iterations(cls, v, values):
        if v:
            return v
        if values.get("spawning") == "independent":
            return 1
        return 30

    @validator("xtc", "pdb", always=True)
    def check_extensions(cls, v, field, values):
        return values.get("traj_name").endswith(f".{field.name}")

    @validator("box_center", always=True)
    def calculate_box_center(cls, v, values):
        if v:
            if ":" in v:
                box_center = [
                    str(coord)
                    for coord in hp.get_coords_from_residue(values.get("system"), v)
                ]
                return "[" + ", ".join(box_center) + "]"
            else:
                return "[" + ",".join([str(coord) for coord in v]) + "]"
        return cls.simulation_params.get("box_center", None)

    @validator("core", always=True)
    def set_core_constraints(cls, v, values):
        if v:
            return v
        return [] if values.get("use_peleffy") else -1

    @validator("pele_exec", always=True)
    def parse_pele_exec(cls, v):
        """
        Sets path to PELE executable following the priority:
            1. Value from YAML file.
            2. Environment variable (PELE_EXEC)
            3. Path_to_PELE + bin/Pele_mpi
        """
        if v:
            return v

        pele_exec = os.environ.get("PELE_EXEC", "")
        default_pele_exec = (
            pele_exec if pele_exec else os.path.join(constants.PELE, "bin/Pele_mpi")
        )
        return default_pele_exec

    @validator("pele_data", always=True)
    def parse_pele_data(cls, v):
        """
        Sets path to PELE Data following the priority:
            1. Value from YAML file.
            2. Environment variable (PELE_DATA)
            3. Path_to_PELE + bin/Pele_mpi
        """
        if v:
            return v

        pele_data = os.environ.get("PELE_DATA", "")
        default_pele_Data = (
            pele_data if pele_data else os.path.join(constants.PELE, "Data")
        )
        return default_pele_Data

    @validator("pele_documents", always=True)
    def parse_pele_documents(cls, v):
        """
        Sets path to PELE Documents following the priority:
            1. Value from YAML file.
            2. Environment variable (PELE_DOCUMENTS)
            3. Path_to_PELE + bin/Pele_mpi
        """
        if v:
            return v

        pele_documents = os.environ.get("PELE_DOCUMENTS", "")
        default_pele_documents = (
            pele_documents
            if pele_documents
            else os.path.join(constants.PELE, "Documents")
        )
        return default_pele_documents

    @validator("license", always=True)
    def parse_pele_license(cls, v):
        """
        Sets path to PELE License following the priority:
            1. Value from YAML file.
            2. Environment variable (PELE_LICENSE)
            3. Path_to_PELE + bin/Pele_mpi
        """
        if v:
            return v

        pele_license = os.environ.get("PELE_LICENSE", "")
        default_pele_license = (
            pele_license if pele_license else os.path.join(constants.PELE, "licenses")
        )
        return default_pele_license

    @validator("singularity_exec", always=True)
    def parse_singularity_executable(cls, v):
        """
        Sets path to Singularity executable following the priority:
            1. Value from YAML file.
            2. Environment variable (SINGULARITY_EXEC)
        """
        return v if v else os.environ.get("SINGULARITY_EXEC", "")

    @validator("water_arg", always=True)
    def parse_water_arg(cls, v, values):
        if values["waters"] == "all_waters":
            return hp.retrieve_all_waters(values["system"])
        return values["waters"]

    @validator("equil_steps")
    def calculate_equilibration_steps(cls, v, values):
        if v:
            return int(v / values["cpus"]) + 1
        return cls.simulation_params.get("equilibration_steps", 1)

    @validator("poses")
    def calculate_poses(cls, v, values):
        if v:
            return v
        return values["cpus"] - 1
