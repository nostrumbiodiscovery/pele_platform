import os
import random
import re
from typing import Any, Optional
from pydantic import BaseModel, validator

from pele_platform.Errors import custom_errors
from pele_platform.Models.utils import Field
from pele_platform.constants import constants


class YamlParserModel(BaseModel):
    class Config:
        allow_population_by_field_name = True
        extra = "allow"

    package: Any = Field()
    adaptive: Any = Field()
    input: Optional[Any] = Field(alias="global_inputs")
    system: Any = Field(default="")
    residue: Any = Field(alias="resname")
    chain: Any = Field()
    hbond: Any = Field(default=[None, None])
    test: Any = Field()
    pele: Any = Field()
    forcefield: Any = Field(
        default="OPLS2005",
        value_from_simulation_params=True,
        simulation_params_default="OPLS2005",  # Redundant, default was already in yaml parser
    )
    verbose: Any = Field()
    anm_freq: Any = Field(
        tests_value=0,
        can_be_falsy=True,
        value_from_simulation_params=True,
        simulation_params_default=4,
    )
    sidechain_freq: Any = Field(
        tests_value=0,
        can_be_falsy=True,
        value_from_simulation_params=True,
        simulation_params_default=2,
    )
    min_freq: Any = Field(
        tests_value=0,
        can_be_falsy=True,
        value_from_simulation_params=True,
        simulation_params_default=1,
    )
    water_freq: Any = Field(
        can_be_falsy=True,
        value_from_simulation_params=True,
        simulation_params_default=1,
    )
    temperature: Any = Field(
        tests_value=10000,
        value_from_simulation_params="temperature",
        simulation_params_default=1500,
    )
    temp: Any = Field(value_from="temperature")
    sidechain_resolution: Any = Field(
        alias="sidechain_res",
        value_from_simulation_params=True,
        simulation_params_default=30,
    )
    steric_trials: Any = Field(
        value_from_simulation_params=True, simulation_params_default=250
    )
    overlap_factor: Any = Field(
        value_from_simulation_params=True, simulation_params_default=0.65
    )
    steering: Any = Field(
        value_from_simulation_params=True, simulation_params_default=0
    )
    solvent: Any = Field(
        value_from_simulation_params=True, simulation_params_default="VDGBNP"
    )
    usesrun: Any = Field()
    spawning: Any = Field(
        value_from_simulation_params="spawning_type",
        simulation_params_default="independent",
    )
    iterations: Any = Field(tests_value=1, value_from_simulation_params=True)
    steps: Any = Field(tests_value=1)
    pele_steps: Any = Field(
        value_from="steps",
        value_from_simulation_params=True,
        simulation_params_default=8,
    )
    cpus: Any = Field(
        tests_value=5, value_from_simulation_params=True, simulation_params_default=60
    )
    density: Any = Field(
        value_from_simulation_params=True,
        simulation_params_default="null",
    )
    cluster_values: Any = Field(
        value_from_simulation_params=True,
        simulation_params_default="[1.75, 2.5, 4, 6]",
    )
    cluster_conditions: Any = Field(
        value_from_simulation_params=True,
        simulation_params_default="[1, 0.6, 0.4, 0.0]",
    )
    simulation_type: Any = Field(
        value_from_simulation_params=True,
        simulation_params_default="pele",
    )
    equilibration: Any = Field()
    clust_type: Any = Field(
        value_from_simulation_params=True,
        simulation_params_default="rmsd",
    )
    eq_steps: Any = Field(alias="equilibration_steps")
    adaptive_restart: Any = Field()
    report_name: Any = Field(alias="report", default="report")
    traj_name: Any = Field(alias="traj", default="trajectory.pdb")
    epsilon: Any = Field(value_from_simulation_params=True, simulation_params_default=0)
    out_in: Any = Field()
    bias_column: Any = Field(
        value_from_simulation_params=True, simulation_params_default=5
    )
    gridres: Any = Field(default=10)
    core: Any = Field(default=-1)
    mtor: Any = Field(alias="maxtorsion", default=4)
    n: int = Field(default=10000)
    template: Any = Field(alias="templates")
    ext_temp: Any = Field(value_from="template")
    rotamers: Any = Field()
    ext_rotamers: Any = Field(value_from="rotamers")
    mae_lig: Any = Field(value_from_simulation_params=True)
    no_ppp: Any = Field(
        alias="skip_preprocess",
        value_from_simulation_params=True,
        simulation_params_default=False,
    )
    skip_prep: Any = Field(
        alias="skip_preprocess",
        value_from_simulation_params=True,
        simulation_params_default=False,
    )
    gaps_ter: Any = Field(
        alias="TERs",
        can_be_falsy=True,
        value_from_simulation_params=True,
        simulation_params_default=False,
    )
    charge_ter: Any = Field(
        alias="charge_ters",
        can_be_falsy=True,
        value_from_simulation_params=True,
        simulation_params_default=False,
    )
    mpi_params: Any = Field()
    nonstandard: Any = Field(
        value_from_simulation_params=True, simulation_params_default=[]
    )
    prepwizard: Any = Field()
    box_radius: Any = Field(value_from_simulation_params=True)
    box_center: Any = Field()
    box: Any = Field()
    native: Any = Field(alias="rmsd_pdb", default=False)
    atom_dist: Any = Field(default_factory=list)
    debug: bool = Field(default=False)
    folder: Any = Field(alias="working_folder")
    output: Any = Field(default="output")
    randomize: Any = Field(
        value_from_simulation_params=True, simulation_params_default=False
    )
    full: Any = Field(alias="global")
    proximityDetection: Any = Field()
    poses: Any = Field()
    precision_glide: Any = Field()
    msm: Any = Field()
    precision: Any = Field()
    clust: Any = Field(alias="exit_clust", tests_value=2)
    restart: Any = Field(
        value_from_simulation_params=True, simulation_params_default="all"
    )
    lagtime: Any = Field()
    msm_clust: Any = Field()
    rescoring: Any = Field()
    in_out: Any = Field()
    in_out_soft: Any = Field()
    exit: Any = Field()
    exit_value: Any = Field()
    exit_condition: Any = Field()
    exit_trajnum: Any = Field()
    waters: Any = Field(value_from_simulation_params=True, simulation_params_default=[])
    water_center: Any = Field()
    water_temp: Any = Field(
        value_from_simulation_params=True, simulation_params_default=5000
    )
    water_overlap: Any = Field(
        value_from_simulation_params=True, simulation_params_default=0.78
    )
    water_constr: Any = Field(
        value_from_simulation_params=True, simulation_params_default=0
    )
    water_trials: Any = Field(
        value_from_simulation_params=True, simulation_params_default=10000
    )
    water_radius: Any = Field(default=6)
    induced_fit_exhaustive: Any = Field()
    induced_fit_fast: Any = Field()
    frag: Any = Field()
    ca_constr: Any = Field(
        can_be_falsy=True,
        value_from_simulation_params=True,
        simulation_params_default=5,
    )
    ca_interval: Any = Field(
        can_be_falsy=True,
        value_from_simulation_params=True,
        simulation_params_default=5,
    )
    one_exit: Any = Field()
    box_type: Any = Field()
    box_metric: Any = Field()
    time: Any = Field()
    nosasa: Any = Field()
    perc_sasa: Any = Field()
    seed: Any = Field(default_factory=lambda: random.randrange(1, 70000))
    pdb: Any = Field()
    log: Any = Field()
    nonrenum: Any = Field()
    pele_exec: Any = Field(default=os.path.join(constants.PELE, "bin/Pele_mpi"))
    pele_data: Any = Field(default=os.path.join(constants.PELE, "Data"))
    pele_documents: Any = Field(default=os.path.join(constants.PELE, "Documents"))
    pca: Any = Field()
    anm_direction: Any = Field(
        value_from_simulation_params=True,
        simulation_params_default="random",
    )
    anm_mix_modes: Any = Field(
        value_from_simulation_params=True,
        simulation_params_default="mixMainModeWithOthersModes",
    )
    anm_picking_mode: Any = Field(
        value_from_simulation_params=True,
        simulation_params_default="RANDOM_MODE",
    )
    anm_displacement: Any = Field(
        value_from_simulation_params=True, simulation_params_default=0.75
    )
    anm_modes_change: Any = Field(
        value_from_simulation_params=True, simulation_params_default=4
    )
    anm_num_of_modes: Any = Field(
        value_from_simulation_params=True, simulation_params_default=6
    )
    anm_relaxation_constr: Any = Field(
        value_from_simulation_params=True, simulation_params_default=0.5
    )
    remove_constraints: Any = Field(
        can_be_falsy=True,
        value_from_simulation_params=True,
        simulation_params_default=False,
    )
    pca_traj: Any = Field()
    perturbation: Any = Field()
    sasa: Any = Field(
        value_from_simulation_params=True, simulation_params_default=constants.SASA
    )
    binding_energy: Any = Field(
        value_from_simulation_params=True, simulation_params_default=constants.BE
    )
    parameters: Any = Field(
        value_from_simulation_params="params", simulation_params_default=True
    )
    analyse: Any = Field(default=True)
    selection_to_perturb: Any = Field(
        value_from_simulation_params=True,
        simulation_params_default=constants.SELECTION_TO_PERTURB,
    )
    mae: Any = Field(default=False)
    constrain_core: Any = Field(value_from_simulation_params=True)
    constrain_core_spring: Any = Field(default=50.0)
    skip_ligand_prep: Any = Field(
        value_from_simulation_params="args.skip_ligand_prep",
        simulation_params_default=[],
    )
    spawning_condition: Any = Field(value_from_simulation_params=True)
    external_constraints: Any = Field(default=[])
    only_analysis: Any = Field(default=False)
    overwrite: Any = Field(alias="overwrite_analysis", default=True)
    analysis_nclust: Any = Field(default=10)
    te_column: Any = Field(default=4)
    be_column: Any = Field(default=5)
    limit_column: Any = Field(default=6)
    com: Any = Field(
        alias="COMligandConstraint",
        value_from_simulation_params="COMligandConstraint",
        simulation_params_default=0,
    )
    pele_license: Any = Field(default=os.path.join(constants.PELE, "licenses"))
    license: Any = Field(value_from="pele_license")
    schrodinger: Any = Field()
    no_check: Any = Field(default=False)
    cleanup: Any = Field(default=False)
    water_empty_selector: Any = Field(default=False)
    polarize_metals: Any = Field(default=False)
    polarization_factor: Any = Field(default=2.0)
    workflow: Any = Field()
    distance: Any = Field()
    permissive_metal_constr: Any = Field(default_factory=list)
    constrain_all_metals: Any = Field(default=False)
    no_metal_constraints: Any = Field(default=False)
    frag_run: Any = Field(default=True)
    frag_core: Any = Field(default=False)
    frag_input: Any = Field(default=False)
    frag_ligands: Any = Field(default=False)
    growing_steps: Any = Field(default=False)
    frag_steps: Any = Field(alias="steps_in_gs", default=False)
    frag_eq_steps: Any = Field(alias="sampling_steps", default=False)
    protocol: Any = Field()
    frag_ai: Any = Field(default=False)
    frag_ai_iterations: Any = Field(default=False)
    chain_core: Any = Field(default=False)
    frag_restart: Any = Field(default=False)
    frag_criteria: Any = Field(default=False)
    frag_output_folder: Any = Field(default=False)
    frag_cluster_folder: Any = Field(default=False)
    frag_library: Any = Field()
    frag_core_atom: Any = Field()
    analysis_to_point: Any = Field()
    # allosteric defined 10, and ppi overwrote it to 25
    n_components: Any = Field(
        tests_value=3,
        value_from_simulation_params=True,
        simulation_params_default=25,
    )
    ppi: Any = Field()
    center_of_interface: Any = Field()
    protein: Any = Field()
    ligand_pdb: Any = Field()
    skip_refinement: Any = Field(default=False)
    n_waters: Any = Field(
        value_from_simulation_params=True, simulation_params_default=0
    )
    allosteric: Any = Field()
    rna: Any = Field()
    gpcr_orth: Any = Field()
    orthosteric_site: Any = Field()
    initial_site: Any = Field()
    final_site: Any = Field()

    @validator("*", pre=True, always=True)
    def set_tests_values(cls, v, values, field):
        if values.get("test"):
            test_value = field.field_info.extra.get("tests_value")
            if test_value is not None:
                return test_value
        return v

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

    @validator("initial_site", "orthosteric_site", "final_site")
    def validate_atom_string(cls, v):
        if v:
            pattern = r"([A-z]\:\d{1,4}\:[A-Z0-9]{1,4})"
            if not re.match(pattern, v):
                raise custom_errors.WrongAtomStringFormat(
                    "Atom string set in {} does not seem to have the right format. It should follow chain:residue "
                    "number:atom name patter, e.g. 'A:105:CA'".format(
                        v
                    )
                )
