import os
import glob
from typing import Any, List
import pele_platform.constants.constants as cs
from pydantic import validator
import pele_platform.Utilities.Helpers.helpers as hp
from pele_platform.Models.utils import Field
from pele_platform.Models.yaml_parser_model import YamlParserModel

LOGFILE = '"simulationLogPath" : "$OUTPUT_PATH/logFile.txt",'


class SimulationParamsModel(YamlParserModel):
    frag_pele: Any = Field()
    complexes: str = Field()
    frag_pele_steps: Any = Field()
    output_path: Any = Field()
    logfile: str = Field(value_from="log", always=True)
    water: str = Field(default="")
    ligand: Any = Field(default=cs.LIGAND)
    spython: Any = Field(default=os.path.join(cs.SCHRODINGER, "utilities/python"))
    lig: Any = Field(value_from="mae_lig")
    sasa_max: Any = Field()
    sasa_min: Any = Field()
    clusters: Any = Field(alias="clust", tests_value=2)
    allow_empty_selectors: Any = Field(value_from="water_empty_selector")
    templates: Any = Field(
        default=os.path.abspath(
            os.path.join(os.path.dirname(os.path.dirname(__file__)), "PeleTemplates")
        )
    )
    xtc: bool = Field()
    pdb: bool = Field()
    constraints: Any = Field()
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
    proximityDetection: bool = Field()

    @validator("*", pre=True, always=True)
    def set_value_from_simulation_parameters(cls, v, field):

        extra = field.field_info.extra
        can_be_falsy = extra.get("can_be_falsy")

        if (can_be_falsy and v is None) or (not can_be_falsy and not v):
            value_from_simulation_params = extra.get("value_from_simulation_params")

            if value_from_simulation_params == True:  # TODO: IS?
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
        if values.get("package") in ["site_finder", "adaptive", "PPI", "gpcr_orth"]:
            return True

    @validator("frag_pele", always=True)
    def set_frag_pele(cls, v, values):
        if values.get("frag_core"):
            return True
        else:
            return False

    @validator("complexes", always=True)
    def set_complexes(cls, v, values):
        return "$PDB" if values.get("frag_pele") else "$COMPLEXES"

    @validator("frag_pele_steps", always=True)
    def set_frag_pele_steps(cls, v, values):
        return "$STEPS" if values.get("frag_pele") else "$PELE_STEPS"

    # @validator("chain")
    # def set_frag_chain(cls, v, values):
    #     if values.get("frag_pele"):
    #         return values.get("frag_core")

    @validator("output_path", always=True)
    def set_output_path(cls, v, values):
        return "$RESULTS_PATH" if values.get("frag_pele") else "$OUTPUT_PATH"

    @validator("input", always=True)
    def set_input_glob(cls, v, values):
        system = values.get("system")
        if system and "*" in system:
            return glob.glob(system)
        return v

    @validator("system", always=True)
    def set_system_glob(cls, v):
        if v and "*" in v:
            return glob.glob(v)[0]
        return v

    @validator("logfile", always=True)
    def set_logfile(cls, v):
        if v:
            return LOGFILE
        return ""

    @validator("system", "residue", "chain", always=True)
    def validate_adaptive_required_fields(cls, v, values, field):
        if values.get("adaptive") and not v:
            raise AssertionError(f"User must define {field.name} to use AdaptivePELE.")
        return v

    @validator("perturbation", always=True)
    def set_perturbation(cls, v, values):
        if v is False or values.get("pca_traj"):
            return ""
        return cls.simulation_params.get("perturbation", cs.PERTURBATION)

    @validator("sidechain_perturbation", always=True)
    def set_sidechain_perturbation(cls, v, values):

        if v or values.get("covalent_residue"):
            return True
        else:
            return False

    @validator("inter_step_logger", always=True)
    def set_interstep_logger(cls, v):
        if v is True:
            return cs.INTERSTEPLOGGER
        return ""

    @validator("minimum_steps", always=True)
    def set_minimum_steps(cls, v):
        if v is True:
            return cs.MINIMUMSTEPS
        return ""

    @validator("proximityDetection", always=True)
    def set_proximityDetection(cls, v):
        if v is False:
            return "false"
        return cls.simulation_params.get("proximityDetection", "true")

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

    @validator("conformation_perturbation", always=True)
    def set_conformation_perturbation(cls, v, values):
        if values.get("ligand_conformations"):
            return cs.CONFORMATION_PERTURBATION
        return ""

    @validator("conformation_freq", always=True)
    def only_with_conformation_perturbation(cls, v, values):
        if values.get("ligand_conformations"):
            return cs.CONFORMATION_FREQUENCY.format(v)
        return ""

    @validator("spython", always=True)
    def check_spython_path(cls, v):
        if not os.path.exists(v):
            return os.path.join(cs.SCHRODINGER, "run")
        return v

    @validator("lig", always=True)
    def set_lig(cls, v, values):
        if v:
            return v
        return "{}.mae".format(values.get("residue"))

    @validator("spawning_condition", always=True)
    def format_spawning_condition(cls, v):
        if v:
            return f'"condition": "{v}",'
        return ""

    @validator("allow_empty_selectors", always=True)
    def format_allow_empty_selectors(cls, v):
        if v:
            '"allowEmptyWaterSelectors": true,'
        return ""

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

    @validator("verbose", "equilibration", always=True)
    def format_verbose(cls, v, field):
        return "true" if v else cls.simulation_params.get(field.name, "false")

    @validator("core", always=True)
    def set_core_constraints(cls, v, values):
        if v:
            return v
        else:
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
        else:
            pele_exec = os.environ.get("PELE_EXEC", "")

            default_pele_exec = (
                pele_exec
                if pele_exec
                else os.path.join(cs.PELE, "bin/Pele_mpi")
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
        else:
            pele_data = os.environ.get("PELE_DATA", "")
            default_pele_Data = (
                pele_data
                if pele_data
                else os.path.join(cs.PELE, "Data")
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
        else:
            pele_documents = os.environ.get("PELE_DOCUMENTS", "")
            default_pele_documents = (
                pele_documents
                if pele_documents
                else os.path.join(cs.PELE, "Documents")
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
        else:
            pele_license = os.environ.get("PELE_LICENSE", "")
            default_pele_license = (
                pele_license
                if pele_license
                else os.path.join(cs.PELE, "licenses")
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
