import os
import random
import re
from typing import Any, Optional, List, Union
from pydantic import BaseModel, validator

from pele_platform.Errors import custom_errors
from pele_platform.Models.utils import Field
from pele_platform.constants import constants


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
    system: str = Field(default="", categories=["General settings"])
    residue: str = Field(alias="resname", categories=["General settings"])
    chain: str = Field(categories=["General settings"])
    hbond: Any = Field(
        default=[None, None], candidate_for_deprecation=True
    )  # Kill it with fire!
    test: bool = Field(categories=["General settings"])

    pele: Any = Field(candidate_for_deprecation=True)
    forcefield: str = Field(
        default="OPLS2005",
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
    sidechain_resolution: int = Field(  # int add validator? 360 %% v == 0
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
    cluster_values: List[
        float
    ] = Field(  # WTF - List[float] but then use validator to make it into a str(v), what about bool?
        value_from_simulation_params=True,
        simulation_params_default="[1.75, 2.5, 4, 6]",
        categories=["Simulation parameters"],
    )
    cluster_conditions: List[float] = Field(
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
    )  # alias="ligand_resolution" ? validator for 360 modulus
    core: int = Field(default=-1, categories=["Ligand preparation"])

    ################################################################################################ start from here

    mtor: int = Field(alias="maxtorsion", default=4, categories=["Ligand preparation"])
    n: int = Field(
        default=10000,
        categories=["Ligand preparation"],
        description="Maximum number of flexible sidechains in the ligand.",
    )
    template: List[str] = Field(
        alias="templates"
    )  # WTF with templates - template - ext_temp???
    ext_temp: List[str] = Field(
        value_from="template", categories=["Ligand preparation"]
    )
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
    )
    charge_ter: bool = Field(
        alias="charge_ters",
        can_be_falsy=True,
        value_from_simulation_params=True,
        simulation_params_default=False,
    )
    mpi_params: Any = Field()  # WTF - what type is that supposed to be?
    nonstandard: List[str] = Field(
        value_from_simulation_params=True, simulation_params_default=[]
    )
    prepwizard: bool = Field()
    box_radius: float = Field(  # should be a float but tests fail --> FLOAT YES
        value_from_simulation_params=True, categories=["Box settings"]
    )
    box_center: Union[List[float], str] = Field(categories=["Box settings"])
    box: Any = Field(categories=["Box settings"])
    native: str = Field(alias="rmsd_pdb", default=False, categories=["Metrics"])
    atom_dist: Union[List[str], List[int]] = Field(
        default_factory=list, categories=["Metrics"]
    )  # can be atom numbers or strings
    debug: bool = Field(default=False, categories=["General settings"])
    folder: str = Field(alias="working_folder", categories=["General settings"])
    output: str = Field(default="output", categories=["General settings"])
    randomize: bool = Field(
        value_from_simulation_params=True,
        simulation_params_default=False,
        categories=["Global Exploration", "Simulation parameters"],
    )
    full: bool = Field(alias="global", categories=["Global Exploration"])
    proximityDetection: bool = Field()  # WTF?
    poses: int = Field()
    precision_glide: Any = Field(candidate_for_deprecation=True)
    msm: Any = Field()
    precision: Any = Field(candidate_for_deprecation=True)  # WTF
    clust: Any = Field(
        alias="exit_clust", tests_value=2, candidate_for_deprecation=True
    )

    # WTF?
    restart: Any = Field(
        value_from_simulation_params=True, simulation_params_default="all"
    )
    lagtime: Any = Field(candidate_for_deprecation=True)  # WTF?
    msm_clust: Any = Field()
    rescoring: bool = Field()
    in_out: bool = Field()
    in_out_soft: bool = Field()
    exit: Any = Field()  # WTF?
    exit_value: float = Field()
    exit_condition: str = Field()
    exit_trajnum: int = Field()
    waters: Union[str, List[str]] = Field(
        value_from_simulation_params=True,
        simulation_params_default=[],
        categories=["Water"],
    )
    water_center: Union[List[float], str] = Field(categories=["Water"])
    water_temp: Any = Field(
        value_from_simulation_params=True,
        simulation_params_default=5000,
        categories=["Water"],
    )
    water_overlap: float = Field(
        value_from_simulation_params=True,
        simulation_params_default=0.78,
        categories=["Water"],
    )
    water_constr: Any = Field(
        value_from_simulation_params=True,
        simulation_params_default=0,
        categories=["Water"],
    )
    water_trials: Any = Field(
        value_from_simulation_params=True,
        simulation_params_default=10000,
        categories=["Water"],
    )
    water_radius: int = Field(
        default=6, categories=["Water"]
    )  # WTF tests expects int but should be a float perhaps
    induced_fit_exhaustive: bool = Field(categories=["Induced fit"])
    induced_fit_fast: bool = Field(categories=["Induced fit"])
    frag: Any = Field(categories=["FragPELE"])
    ca_constr: int = Field(
        can_be_falsy=True,
        value_from_simulation_params=True,
        simulation_params_default=5,
        categories=["Constraints"],
    )
    ca_interval: float = Field(
        can_be_falsy=True,
        value_from_simulation_params=True,
        simulation_params_default=5,
        categories=["Constraints"],
    )
    one_exit: Any = Field()
    box_type: Any = Field()
    box_metric: Any = Field()
    time: Any = Field()
    nosasa: Any = Field()
    perc_sasa: Any = Field()
    seed: int = Field(default_factory=lambda: random.randrange(1, 70000))
    pdb: Any = Field()
    log: Any = Field()
    nonrenum: Any = Field()
    pele_exec: str = Field(default=os.path.join(constants.PELE, "bin/Pele_mpi"))
    pele_data: str = Field(default=os.path.join(constants.PELE, "Data"))
    pele_documents: str = Field(default=os.path.join(constants.PELE, "Documents"))
    pca: Any = Field()
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
    anm_modes_change: Any = Field(
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
    pca_traj: List[str] = Field()
    perturbation: Any = Field()
    sasa: str = Field(
        value_from_simulation_params=True, simulation_params_default=constants.SASA
    )
    binding_energy: str = Field(
        value_from_simulation_params=True, simulation_params_default=constants.BE
    )
    parameters: str = Field(
        value_from_simulation_params="params", simulation_params_default=True
    )
    analyse: Any = Field(default=True)
    selection_to_perturb: str = Field(
        value_from_simulation_params=True,
        simulation_params_default=constants.SELECTION_TO_PERTURB,
    )
    mae: Any = Field(default=False)
    constrain_core: Any = Field(value_from_simulation_params=True)
    constrain_core_spring: int = Field(default=50.0, categories=["Constraints"])
    skip_ligand_prep: List[str] = Field(
        value_from_simulation_params="args.skip_ligand_prep",
        simulation_params_default=[],
        categories=["Ligand preparation"],
    )
    spawning_condition: Any = Field(value_from_simulation_params=True)
    external_constraints: List[str] = Field(default=[], categories=["Constraints"])
    only_analysis: bool = Field(default=False, categories=["Analysis"])
    overwrite: bool = Field(alias="overwrite_analysis", default=True)
    analysis_nclust: int = Field(default=10)
    te_column: int = Field(default=4)
    be_column: int = Field(default=5)
    limit_column: int = Field(default=6)
    com: Any = Field(
        alias="COMligandConstraint",
        value_from_simulation_params="COMligandConstraint",
        simulation_params_default=0,
    )
    pele_license: str = Field(default=os.path.join(constants.PELE, "licenses"))
    license: Any = Field(value_from="pele_license")
    schrodinger: str = Field()
    no_check: bool = Field(default=False)
    cleanup: bool = Field(default=False, categories=["FragPELE"], description="Automatically cleans up fragment files, only applicable to FragPELE.")
    water_empty_selector: Any = Field(default=False, categories=["Water"])
    polarize_metals: bool = Field(default=False, categories=["Metals"])
    polarization_factor: float = Field(default=2.0, categories=["Metals"])
    workflow: List[Any] = Field(categories=["Custom workflows"])
    distance: float = Field(categories=["Custom workflows"])
    permissive_metal_constr: Any = Field(
        default_factory=list, categories=["Constraints"]
    )
    constrain_all_metals: Any = Field(default=False, categories=["Constraints"])
    no_metal_constraints: Any = Field(default=False, categories=["Constraints"])
    frag_run: bool = Field(default=True, categories=["FragPELE"])
    frag_core: str = Field(categories=["FragPELE"])
    frag_input: str = Field(default=False, categories=["FragPELE"])
    frag_ligands: Any = Field(default=False, categories=["FragPELE"])
    growing_steps: int = Field(default=False, categories=["FragPELE"])
    frag_steps: int = Field(alias="steps_in_gs", default=False, categories=["FragPELE"])
    frag_eq_steps: int = Field(
        alias="sampling_steps", default=False, categories=["FragPELE"]
    )
    protocol: str = Field(categories=["FragPELE"])
    frag_ai: bool = Field(default=False, categories=["FragPELE"])
    frag_ai_iterations: int = Field(default=False, categories=["FragPELE"])
    chain_core: str = Field(default=False, categories=["FragPELE"])
    frag_restart: bool = Field(default=False, categories=["FragPELE"])
    frag_criteria: str = Field(default=False, categories=["FragPELE"])
    frag_output_folder: str = Field(default=False, categories=["FragPELE"])
    frag_cluster_folder: str = Field(default=False, categories=["FragPELE"])
    frag_library: str = Field(categories=["FragPELE"])
    frag_core_atom: str = Field(
        categories=["FragPELE"]
    )  # we could add a validator as well, e.g. C3-H2
    analysis_to_point: Optional[List[float]] = Field(categories=["FragPELE"])

    n_components: int = Field(
        tests_value=3,
        value_from_simulation_params=True,
        simulation_params_default=10,
        categories=["PPI", "Allosteric"],
    )
    ppi: bool = Field(categories=["PPI"])
    center_of_interface: str = Field(categories=["PPI"])
    protein: str = Field(categories=["PPI"])
    ligand_pdb: str = Field(categories=["PPI"])
    skip_refinement: bool = Field(default=False, categories=["PPI", "Allosteric"])
    n_waters: int = Field(
        value_from_simulation_params=True,
        simulation_params_default=0,
        categories=["Water"],
    )
    allosteric: bool = Field(categories=["Allosteric"])
    rna: bool = Field(categories=["RNA"])
    gpcr_orth: bool = Field(categories=["GPCR"])
    orthosteric_site: str = Field(categories=["GPCR"])
    initial_site: str = Field(categories=["GPCR", "Out in"])
    final_site: str = Field(categories=["GPCR"])

    @validator("*", pre=True, always=True)
    def set_tests_values(cls, v, values, field):
        if values.get("test"):
            test_value = field.field_info.extra.get("tests_value")
            if test_value is not None:
                return test_value
        return v

    # TODO: do we automatically swap floats to ints?

    @validator("*", pre=True, always=True)
    def set_value_from(cls, v, values, field):
        value_from = field.field_info.extra.get("value_from")
        if value_from is not None:
            return values.get(value_from)
        return v

    @validator("system", "mae_lig")
    def construct_path(cls, v):
        return os.path.abspath(v) if v else None

    @validator("residue")
    def validate_residue(cls, v):
        if v == "UNK":
            raise custom_errors.LigandNameNotSupported(
                "'UNK' ligand name is not supported, please rename it, e.g. 'LIG'."
            )
        return v

    @validator("usesrun", always=True)
    def set_usesrun(cls, v):
        return bool(os.environ.get("SRUN", v))

    @validator("initial_site", "orthosteric_site", "final_site", "center_of_interface")
    def validate_atom_string(cls, v):
        if v:
            pattern = r"([A-z]\:\d{1,4}\:[A-Z0-9]{1,4})"
            if not re.match(pattern, v):
                raise custom_errors.WrongAtomStringFormat(
                    "Atom string set in {} does not seem to have the right format. It should follow chain:residue "
                    "number:atom name patter, e.g. 'A:105:CA'".format(v)
                )
        return v
