import os
import random
import re
from typing import Any, Optional, List, Union
from pydantic import BaseModel, validator

from pele_platform.Errors import custom_errors
from pele_platform.Models.utils import Field
from pele_platform.constants import constants


def generate_random_seed():
    return random.randrange(1, 70000)


class YamlParserModel(BaseModel):
    class Config:
        allow_population_by_field_name = True
        extra = "allow"
        # validate_all TODO: consider validating all defaults as well

    package: str = Field(categories=["General settings"])
    adaptive: str = Field(categories=["General settings"])
    input: Optional[Any] = Field(
        alias="global_inputs", categories=["General settings"]
    )  # check what this is
    system: Optional[str] = Field(default="", categories=["General settings"])
    residue: str = Field(alias="resname", categories=["General settings"])
    chain: Optional[str] = Field(categories=["General settings"])
    hbond: Any = Field(
        default=[None, None], candidate_for_deprecation=True
    )  # Kill it with fire!
    test: bool = Field(categories=["General settings"])

    pele: Any = Field(candidate_for_deprecation=True)
    forcefield: str = Field(

        value_from_simulation_params=True,
        simulation_params_default="OPLS2005",  # Redundant, default was already in yaml parser
        categories=["Simulation parameters"],
    )
    verbose: bool = Field(categories=["General settings"])
    anm_freq: int = Field(
        tests_value=0,
        can_be_falsy=True,
        value_from_simulation_params=True,
        simulation_params_default=4,
        categories=["Simulation parameters"],
    )
    sidechain_freq: int = Field(
        tests_value=0,
        can_be_falsy=True,
        value_from_simulation_params=True,
        simulation_params_default=2,
        categories=["Simulation parameters"],
    )
    min_freq: int = Field(
        tests_value=0,
        can_be_falsy=True,
        value_from_simulation_params=True,
        simulation_params_default=1,
        categories=["Simulation parameters"],
    )
    water_freq: int = Field(
        can_be_falsy=True,
        value_from_simulation_params=True,
        simulation_params_default=1,
        categories=["Simulation parameters"],
    )
    temperature: int = Field(  # check if PELE takes float
        tests_value=10000,
        value_from_simulation_params="temperature",
        simulation_params_default=1500,
        categories=["Simulation parameters"],
    )
    temp: int = Field(
        value_from="temperature", candidate_for_deprecation=True
    )  # probably get rid of it
    sidechain_resolution: int = Field(
        alias="sidechain_res",
        value_from_simulation_params=True,
        simulation_params_default=30,
        categories=["Simulation parameters"],
    )
    steric_trials: int = Field(
        value_from_simulation_params=True,
        simulation_params_default=250,
        categories=["Simulation parameters"],
    )
    overlap_factor: float = Field(  # validator between 0 and 1
        value_from_simulation_params=True,
        simulation_params_default=0.65,
        categories=["Simulation parameters"],
    )
    steering: int = Field(  # WTF type? not sure what that is, int based on test inputs
        value_from_simulation_params=True,
        simulation_params_default=0,
        categories=["Simulation parameters"],
        description="Number of translations in the same direction.",  # check PELE++ docs
    )
    solvent: str = Field(
        value_from_simulation_params=True,
        simulation_params_default="VDGBNP",
        categories=["Simulation parameters"],
    )
    usesrun: bool = Field(categories=["Miscellaneous"])
    spawning: str = Field(  # validator for available options
        value_from_simulation_params="spawning_type",
        simulation_params_default="independent",
        categories=["Simulation parameters"],
    )

    iterations: int = Field(
        tests_value=1,
        value_from_simulation_params=True,
        categories=["Simulation parameters"],
    )
    steps: int = Field(
        tests_value=1,
        categories=["Simulation parameters"],
        candidate_for_deprecation=True,
    )
    pele_steps: int = Field(
        value_from="steps",
        value_from_simulation_params=True,
        simulation_params_default=8,
        categories=["Simulation parameters"],
    )
    cpus: int = Field(
        tests_value=5,
        value_from_simulation_params=True,
        simulation_params_default=60,
        categories=["Simulation parameters"],
    )
    density: str = Field(  # validator to make sure user chooses from available options
        value_from_simulation_params=True,
        simulation_params_default="null",
        categories=["Simulation parameters"],
    )
    cluster_values: Union[
        List[float], str
    ] = Field(  # WTF - List[float] but then use validator to make it into a str(v), what about bool?
        value_from_simulation_params=True,
        simulation_params_default="[1.75, 2.5, 4, 6]",
        categories=["Simulation parameters"],
    )
    cluster_conditions: Union[List[float], str] = Field(
        value_from_simulation_params=True,
        simulation_params_default="[1, 0.6, 0.4, 0.0]",
        categories=["Simulation parameters"],
    )
    simulation_type: str = (
        Field(  # AdaptivePELE -> pele or md, deafult to pele (correct)
            value_from_simulation_params=True,
            simulation_params_default="pele",
            categories=["Simulation parameters"],
        )
    )
    equilibration: bool = Field(categories=["Simulation parameters"])
    clust_type: str = Field(  # validate if the user chose one of the available options, check Adaptive docs
        value_from_simulation_params=True,
        simulation_params_default="rmsd",
        categories=["Simulation parameters"],
    )
    eq_steps: int = Field(
        alias="equilibration_steps", categories=["Simulation parameters"]
    )
    adaptive_restart: bool = Field()
    report_name: str = Field(
        alias="report",
        default="report",
        categories=["Simulation parameters"],
        candidate_for_deprecation=True,
    )
    traj_name: str = Field(  # for XTC (validator? os.path.splitext -- > XTC or PDB)
        alias="traj", default="trajectory.pdb", categories=["Simulation parameters"]
    )

    epsilon: float = Field(
        value_from_simulation_params=True,
        simulation_params_default=0,
        categories=["Simulation parameters"],
    )
    out_in: bool = Field(categories=["Out in"])
    bias_column: int = Field(
        value_from_simulation_params=True,
        simulation_params_default=5,
        categories=["Simulation parameters"],
    )
    gridres: int = Field(
        default=10, categories=["Ligand preparation"]
    )  # alias="ligand_resolution"
    use_peleffy: bool = Field(
        categories=["Ligand preparation"], can_be_falsy=True, default=False
    )
    core: Union[int, List[str]] = Field(categories=["Ligand preparation"])

    ################################################################################################ start from here

    mtor: int = Field(alias="maxtorsion", default=4, categories=["Ligand preparation"])
    n: int = Field(
        default=10000,
        categories=["Ligand preparation"],
        description="Maximum number of flexible side chains in the ligand.",
    )

    templates: List[str] = Field(
        categories=["Ligand preparation"],
        description="List of template file containing ligand forcefield parameters.",
    )

    solvent_template: str = Field()

    rotamers: List[str] = Field(categories=["Ligand preparation"])
    ext_rotamers: List[str] = Field(value_from="rotamers")
    mae_lig: str = Field(
        value_from_simulation_params=True,
        categories=["Ligand preparation"],
        description="Maestro file to extract quantum charges.",
    )
    no_ppp: bool = Field(  # WTF skip_preprocess, no_ppp, skip_prep
        alias="skip_preprocess",
        value_from_simulation_params=True,
        simulation_params_default=False,
        categories=["Protein preparation"],
        candidate_for_deprecation=True,
    )
    skip_prep: bool = Field(
        alias="skip_preprocess",
        value_from_simulation_params=True,
        simulation_params_default=False,
        categories=["Protein preparation"],
    )
    gaps_ter: bool = Field(
        alias="TERs",
        can_be_falsy=True,
        value_from_simulation_params=True,
        simulation_params_default=False,
        categories=["Protein preparation"],
    )
    charge_ter: bool = Field(
        alias="charge_ters",
        can_be_falsy=True,
        value_from_simulation_params=True,
        simulation_params_default=False,
        categories=["Protein preparation"],
    )
    mpi_params: str = Field(categories=["General settings"])
    nonstandard: List[str] = Field(
        value_from_simulation_params=True,
        simulation_params_default=[],
        categories=["Protein preparation"],
    )
    prepwizard: bool = Field()
    box_radius: float = Field(  # should be a float but tests fail --> FLOAT YES
        value_from_simulation_params=True, categories=["Box settings"]
    )
    box_center: Union[List[float], str] = Field(categories=["Box settings"])
    box: Any = Field(categories=["Box settings"])
    native: str = Field(alias="rmsd_pdb", default="", categories=["Metrics"])
    atom_dist: List[str] = Field(
        default_factory=list, categories=["Metrics"]
    )  # deprecate atom numbers
    debug: bool = Field(default=False, categories=["General settings"])
    folder: str = Field(alias="working_folder", )
    output: str = Field(default="output", categories=["General settings"])
    randomize: bool = Field(
        value_from_simulation_params=True,
        simulation_params_default=False,
        categories=["Global Exploration", "Simulation parameters"],
    )
    full: bool = Field(alias="global", categories=["Global Exploration"])
    proximityDetection: bool = Field()  # description from PELE++
    poses: int = Field()
    msm: bool = Field(candidate_for_deprecation=True, categories=["MSM"])
    clust: int = Field(
        alias="exit_clust", tests_value=2, candidate_for_deprecation=True
    )

    restart: bool = Field(categories=["General settings"])
    lagtime: Any = Field(candidate_for_deprecation=True, categories=["MSM"])  # WTF?
    msm_clust: Any = Field(candidate_for_deprecation=True, categories=["MSM"])
    rescoring: bool = Field(categories=["Rescoring"])
    in_out: bool = Field(categories=["In Out"])
    in_out_soft: bool = Field(categories=["In Out"])
    exit: bool = Field(categories=["In Out"])
    exit_value: float = Field(categories=["In Out"])
    exit_condition: str = Field(categories=["In Out"])
    exit_trajnum: int = Field(categories=["In Out"])
    waters: Union[str, List[str]] = Field(
        value_from_simulation_params=True,
        simulation_params_default=[],
        categories=["Water"],
    )
    water_center: Union[List[float], str] = Field(categories=["Water"])
    water_temp: float = Field(
        value_from_simulation_params=True,
        simulation_params_default=5000,
        categories=["Water"],
    )
    water_overlap: float = Field(
        value_from_simulation_params=True,
        simulation_params_default=0.78,
        categories=["Water"],
    )
    water_constr: float = Field(
        value_from_simulation_params=True,
        simulation_params_default=0,
        categories=["Water"],
    )
    water_trials: int = Field(
        value_from_simulation_params=True,
        simulation_params_default=10000,
        categories=["Water"],
    )
    water_radius: float = Field(
        default=6.0, categories=["Water"]
    )
    induced_fit_exhaustive: bool = Field(categories=["Induced fit"])
    induced_fit_fast: bool = Field(categories=["Induced fit"])
    frag: bool = Field(categories=["FragPELE"], candidate_for_deprecation=True)
    ca_constr: float = Field(
        value_from_simulation_params=True,
        simulation_params_default=0.5,
        categories=["Constraints"],
    )
    ca_interval: int = Field(
        value_from_simulation_params=True,
        simulation_params_default=10,
        categories=["Constraints"],
    )
    one_exit: Any = Field(candidate_for_deprecation=True)
    box_type: str = Field(categories=["Box settings"])
    box_metric: Any = Field(candidate_for_deprecation=True)
    time: Any = Field(candidate_for_deprecation=True)
    nosasa: Any = Field(candidate_for_deprecation=True)
    perc_sasa: Any = Field(candidate_for_deprecation=True)
    seed: int = Field(default_factory=generate_random_seed, tests_value=12345)
    pdb: str = Field(categories=["RNA"])
    log: Any = Field(candidate_for_deprecation=True)
    nonrenum: Any = Field(candidate_for_deprecation=True)
    pele_exec: str = Field()
    pele_data: str = Field()
    pele_documents: str = Field(default=os.path.join(constants.PELE, "Documents"))
    pca: str = Field(description="Path to PCA file", categories=["PCA"])
    anm_direction: str = Field(
        value_from_simulation_params=True,
        simulation_params_default="random",
    )
    anm_mix_modes: str = Field(
        value_from_simulation_params=True,
        simulation_params_default="mixMainModeWithOthersModes",
    )
    anm_picking_mode: str = Field(
        value_from_simulation_params=True,
        simulation_params_default="RANDOM_MODE",
    )
    anm_displacement: float = Field(
        value_from_simulation_params=True, simulation_params_default=0.75
    )
    anm_modes_change: int = Field(
        value_from_simulation_params=True, simulation_params_default=4
    )
    anm_num_of_modes: int = Field(
        value_from_simulation_params=True, simulation_params_default=6
    )
    anm_relaxation_constr: float = Field(
        value_from_simulation_params=True, simulation_params_default=0.5
    )
    remove_constraints: bool = Field(
        can_be_falsy=True,
        value_from_simulation_params=True,
        simulation_params_default=False,
    )
    pca_traj: Union[List[str], str] = Field()
    perturbation: bool = Field(categories=["Advanced"])
    sasa: str = Field(
        value_from_simulation_params=True,
        simulation_params_default=constants.SASA,
        categories=["Advanced"],
    )
    binding_energy: str = Field(
        value_from_simulation_params=True,
        simulation_params_default=constants.BE,
        categories=["Advanced"],
    )
    parameters: str = Field(
        value_from_simulation_params="params",
        simulation_params_default=True,
        categories=["Advanced"],
    )
    analyse: bool = Field(default=True)
    selection_to_perturb: str = Field(
        value_from_simulation_params=True,
        simulation_params_default=constants.SELECTION_TO_PERTURB,
        categories=["Advanced"],
    )
    mae: bool = Field(default=False)
    constrain_core: str = Field(
        value_from_simulation_params=True,
        description="String of SMILES or SMARTS to constrain.",
        categories=["Constraints"],
    )
    constrain_core_spring: float = Field(default=50.0, categories=["Constraints"])
    skip_ligand_prep: List[str] = Field(
        value_from_simulation_params="args.skip_ligand_prep",
        simulation_params_default=[],
        categories=["Ligand preparation"],
    )
    spawning_condition: str = Field(
        value_from_simulation_params=True, description="For min or maximising epsilon"
    )
    external_constraints: List[str] = Field(default=[], categories=["Constraints"])
    only_analysis: bool = Field(default=False, categories=["Analysis"])
    overwrite: bool = Field(alias="overwrite_analysis", default=True)
    analysis_nclust: int = Field(default=10)
    te_column: int = Field(default=4)
    be_column: int = Field(default=5)
    limit_column: int = Field(default=6)
    com: float = Field(
        alias="COMligandConstraint",
        value_from_simulation_params="COMligandConstraint",
        simulation_params_default=0,
        categories=["FragPELE"],
    )
    pele_license: str = Field(default=os.path.join(constants.PELE, "licenses"))
    license: str = Field(value_from="pele_license")
    schrodinger: str = Field()
    no_check: bool = Field(default=False)
    cleanup: bool = Field(
        default=False,
        categories=["FragPELE"],
        description="Automatically cleans up fragment files, only applicable to FragPELE.",
    )
    water_empty_selector: bool = Field(default=False, categories=["Water"])
    polarize_metals: bool = Field(default=False, categories=["Metals"])
    polarization_factor: float = Field(default=2.0, categories=["Metals"])
    workflow: List[Any] = Field(categories=["Custom workflows"])
    distance: float = Field(categories=["Custom workflows"])
    permissive_metal_constr: bool = Field(categories=["Constraints"])
    constrain_all_metals: bool = Field(default=False, categories=["Constraints"])
    no_metal_constraints: bool = Field(default=False, categories=["Constraints"])
    frag_run: bool = Field(default=True, categories=["FragPELE"])
    frag_core: str = Field(categories=["FragPELE"])
    frag_input: str = Field(default=None, categories=["FragPELE"])
    frag_ligands: Any = Field(default=None, categories=["FragPELE"])
    growing_steps: int = Field(default=6, categories=["FragPELE"])
    frag_steps: int = Field(alias="steps_in_gs", default=3, categories=["FragPELE"])
    frag_eq_steps: int = Field(
        alias="sampling_steps", default=20, categories=["FragPELE"]
    )
    protocol: str = Field(default="", categories=["FragPELE"])
    frag_ai: bool = Field(default=False, categories=["FragPELE"])
    frag_ai_iterations: int = Field(default=False, categories=["FragPELE"])
    chain_core: str = Field(default="L", categories=["FragPELE"])
    frag_restart: str = Field(default="", categories=["FragPELE"])
    frag_criteria: str = Field(default="Binding Energy", categories=["FragPELE"])
    frag_output_folder: str = Field(default="growing_steps", categories=["FragPELE"])
    frag_cluster_folder: str = Field(default="clustering_PDBs", categories=["FragPELE"])
    frag_library: str = Field(default=None, categories=["FragPELE"])
    frag_core_atom: Union[str, None] = Field(categories=["FragPELE"], default=None)
    analysis_to_point: Optional[List[float]] = Field(categories=["FragPELE"])
    fragment_atom: str = Field(default=None, categories=["FragPELE"])
    frag_restart_libraries: bool = Field(default=None, categories=["FragPELE"])

    n_components: int = Field(
        tests_value=3,
        value_from_simulation_params=True,
        simulation_params_default=10,
        categories=["PPI", "Site finder"],
        candidate_for_deprecation=True,
    )
    ppi: bool = Field(categories=["PPI"])
    center_of_interface: str = Field(categories=["PPI"])
    protein: str = Field(categories=["PPI"])
    ligand_pdb: str = Field(categories=["PPI"])
    skip_refinement: bool = Field(default=False, categories=["PPI", "Site finder"])
    n_waters: int = Field(
        value_from_simulation_params=True,
        simulation_params_default=0,
        categories=["Water"],
    )
    site_finder: bool = Field(categories=["Site finder"])
    rna: bool = Field(categories=["RNA"])
    gpcr_orth: bool = Field(categories=["GPCR"])
    orthosteric_site: str = Field(categories=["GPCR"])
    initial_site: str = Field(categories=["GPCR", "Out in"])
    final_site: str = Field(categories=["GPCR"])

    top_clusters_criterion: str = Field(
        categories=["Analysis"], default="interaction_25_percentile"
    )
    interaction_restrictions: List[dict] = Field(
        categories=["Interaction restrictions"]
    )

    charge_parametrization_method: str = Field(categories=["Ligand preparation"])

    exclude_terminal_rotamers: bool = Field(
        default=True, categories=["Ligand preparation"]
    )

    singularity_exec: str = Field(categories=["General settings"])

    terminal_constr: float = Field(
        value_from_simulation_params="terminal_constr", simulation_params_default=5.0
    )

    covalent_residue: str = Field(categories=["Covalent docking"])

    nonbonding_radius: float = Field(categories=["Covalent docking"], value_from_simulation_params=True,
                                     simulation_params_default=20.0)

    perturbation_trials: int = Field(categories=["Covalent docking"], value_from_simulation_params=True,
                                     simulation_params_default=10)

    refinement_angle: float = Field(categories=["Covalent docking"], value_from_simulation_params=True,
                                    simulation_params_default=10.0)

    covalent_docking_refinement: bool = Field(categories=["Covalent docking"])

    ligand_conformations: str = Field(categories=["Ligand conformations"])

    conformation_freq: int = Field(
        categories=["Ligand conformations"],
        value_from_simulation_params=True,
        default=4,
    )

    overlap_factor_conformation: float = Field(categories=["Ligand conformations"])

    inter_step_logger: bool = Field(categories=["General settings"])

    minimum_steps: bool = Field()

    site_finder_global: bool = Field()

    site_finder_local: bool = Field()

    # TODO: Add validators for all analysis flags
    kde: bool = Field(categories=["Analysis"])
    kde_structs: int = Field(categories=["Analysis"], default=1000)
    plot_filtering_threshold: float = Field(categories=["Analysis"])
    clustering_filtering_threshold: float = Field(
        categories=["Analysis"], default=0.25, can_be_falsy=True
    )
    clustering_method: str = Field(
        categories=["Analysis"],
        value_from_simulation_params=True,
        simulation_params_default="meanshift",
    )
    cluster_representatives_criterion: str = Field(
        categories=["Analysis"],
        value_from_simulation_params=True,
        simulation_params_default="interaction_5_percentile",
    )
    bandwidth: float = Field(
        categories=["Analysis"],
        value_from_simulation_params=True,
        simulation_params_default=2.5,
    )
    max_top_clusters: int = Field(categories=["Analysis"], can_be_falsy=True, default=8)
    min_population: float = Field(
        categories=["Analysis"],
        value_from_simulation_params=True,
        simulation_params_default=0.01,
    )
    max_top_poses: int = Field(categories=["Analysis"], can_be_falsy=True, default=100)

    saturated_mutagenesis: bool = Field(categories=["Saturated mutagenesis"])
    cpus_per_mutation: int = Field(
        categories=["Saturated mutagenesis"], tests_value=2
    )  # TODO: needs validator
    constraint_level: int = Field(categories=["Constraints"])

    @validator("*", pre=True, always=True)
    def set_tests_values(cls, v, values, field):
        """
        Adjusts values of parameters when running tests.
        """
        if values.get("test"):
            test_value = field.field_info.extra.get("tests_value")
            if test_value is not None:
                return test_value
        return v

    @validator("*", pre=True, always=True)
    def set_value_from(cls, v, values, field):
        """
        Gets value from an a duplicated parameter, e.g. steps and pele_steps.
        """
        value_from = field.field_info.extra.get("value_from")
        if value_from is not None:
            return values.get(value_from)
        return v

    @validator("system", "mae_lig")
    def construct_path(cls, v):
        """
        Gets absolute path to the file.
        """
        return os.path.abspath(v) if v else None

    @validator("residue")
    def validate_residue_name(cls, v):
        """
        Checks the residue name for unsupported 'UNK' value.
        """
        if v == "UNK":
            raise custom_errors.LigandNameNotSupported(
                "'UNK' ligand name is not supported, please rename it, e.g. 'LIG'."
            )
        return v

    @validator("usesrun", always=True)
    def set_usesrun(cls, v):
        """
        Gets SRUN from the environment.
        """
        return bool(os.environ.get("SRUN", v))

    @validator(
        "initial_site",
        "orthosteric_site",
        "final_site",
        "center_of_interface",
        "atom_dist",
    )
    def validate_atom_string(cls, v):
        """
        Checks if a list of strings fits a regex pattern describing an atom.
        """
        pattern = r"([A-z]\:\d{1,4}\:[A-Z0-9]{1,4})"

        if v:
            v_list = [v] if not isinstance(v, list) else v

            for string in v_list:
                if not string.isdigit() and not re.match(pattern, string):
                    raise custom_errors.WrongAtomStringFormat(
                        "Atom string set in {} does not seem to have the right format. It should follow chain:residue "
                        "number:atom name pattern, e.g. 'A:105:CA'.".format(string)
                    )
        return v

    @validator("covalent_residue")
    def validated_residue_string(cls, v):
        """
        Checks if the residue string matches a regex patterns. Correct format example: "A:123".
        """
        pattern = r"(^[A-z]\:\d{1,4}$)"

        if v:
            if not v.isdigit() and not re.match(pattern, v):
                raise custom_errors.WrongAtomStringFormat(
                    "Residue string set in {} does not seem to have the right format. "
                    "It should follow chain:residue_number pattern, e.g. 'A:105'.".format(
                        v
                    )
                )
        return v

    @validator("frag_core_atom")
    def validate_frag_core_atom(cls, v):
        """
        Checks if the string fits a regex pattern indicating the core atom,e.g. C3-H2.
        """
        pattern = r"([A-Z]{1,2}\d{1,2}\-H[A-Z]?\d{1,2})"

        if v:
            if not re.match(pattern, v):
                raise custom_errors.WrongAtomStringFormat(
                    f"Atom string set in {v} does not seem to have the right format. It should follow C atom name: H "
                    "atom name format.")
        return v

    @validator("sidechain_resolution", "gridres")
    def check_divisibility(cls, v):
        if v and not 360 % v == 0:
            raise ValueError(
                "The value should be easily multiplied to obtain 360, e.g. 30, 45 or 10 would be valid."
            )
        return v

    @validator("top_clusters_criterion")
    def check_selection_criterion(cls, v):
        if v and v not in constants.metric_top_clusters_criterion.keys():
            raise ValueError(
                f"Selected criterion value {v} is invalid. Please choose one of: {constants.metric_top_clusters_criterion.keys()}"
            )
        return v

    # TODO: Add validator for what's inside interaction restrictions

    @validator("constraint_level")
    def parse_constraint_level(cls, v, values):
        if v:
            if v not in (0, 1, 2, 3):
                raise ValueError(f"Invalid constraint level {v}, should be 0, 1, 2 or 3.")
        return v

    @validator("terminal_constr", "ca_constr")
    def assert_positive_integer(cls, v):
        if v and v < 0:
            raise ValueError(f"Flag {v} should have a positive value.")
        return v
