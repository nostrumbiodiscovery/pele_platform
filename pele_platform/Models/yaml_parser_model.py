import os
import random
import re
from typing import Any, Optional, List, Union
from pydantic import BaseModel, validator

from pele_platform.Errors import custom_errors
from pele_platform.Models.utils import Field
from pele_platform.constants import constants
from pele_platform.Utilities.Helpers import helpers


def generate_random_seed():
    return random.randrange(1, 70000)


class YamlParserModel(BaseModel):
    class Config:
        allow_population_by_field_name = True
        extra = "allow"

    adaptive: bool = Field(
        description="Run AdaptivePELE using custom parameters.", categories=["Advanced"]
    )
    input: List[str] = Field(
        categories=["General settings"],
        description="Paths to the input files in PDB format (as an alternative to 'system' flag).",
    )
    system: Optional[str] = Field(
        default="",
        categories=["General settings"],
        description="Path to the input system in PBD format.",
    )
    residue: str = Field(
        alias="resname",
        categories=["General settings"],
        description="Residue name of the ligand to be perturbed.",
    )
    chain: str = Field(
        categories=["General settings"],
        description="Chain ID of the ligand to be perturbed.",
    )
    test: bool = Field(
        categories=["General settings"],
        description="Activates test mode to check if the simulation runs as expected. Control files "
        "from the test should never be used for production simulations - temperature, ANM "
        "and minimization parameters are set to unreasonable values to make it run faster.",
        can_be_falsy=True,
    )
    external_templates: Union[str, List[str]] = Field(
        alias="templates",
        default_factory=list,
        categories=["Ligand parametrization"],
        description="Paths to custom template files for hetero molecules.",
    )
    external_rotamers: Union[str, List[str]] = Field(
        alias="rotamers",
        default_factory=list,
        categories=["Ligand parametrization"],
        description="Paths to custom rotamer files for hetero molecules.",
    )
    forcefield: str = Field(
        value_from_simulation_params=True,
        simulation_params_default="OPLS2005",
        categories=["Ligand parametrization"],
        description="Forcefield used to parametrize all hetero molecules.",
    )
    verbose: bool = Field(
        value_from_simulation_params=True,
        can_be_falsy=True,
        simulation_params_default=False,
        categories=["General settings"],
        description="Activates verbose mode in PELE.",
    )
    anm_freq: int = Field(
        tests_value=0,
        can_be_falsy=True,
        value_from_simulation_params=True,
        simulation_params_default=4,
        categories=["Simulation parameters"],
        description="Frequency at which ANM is performed.",
    )
    sidechain_freq: int = Field(
        tests_value=0,
        can_be_falsy=True,
        value_from_simulation_params=True,
        simulation_params_default=2,
        categories=["Simulation parameters"],
        description="Frequency at which side chain sampling is performed.",
    )
    min_freq: int = Field(
        tests_value=0,
        can_be_falsy=True,
        value_from_simulation_params=True,
        simulation_params_default=1,
        categories=["Simulation parameters"],
        description="Frequency at which minimization is performed.",
    )
    water_freq: int = Field(
        can_be_falsy=True,
        value_from_simulation_params=True,
        simulation_params_default=1,
        categories=["Simulation parameters"],
        description="Frequency at which water perturbation is performed.",
    )
    temperature: int = Field(  # check if PELE takes float
        tests_value=10000,
        value_from_simulation_params="temperature",
        simulation_params_default=1500,
        categories=["Simulation parameters"],
        description="Temperature of the simulation.",
    )
    sidechain_resolution: int = Field(
        alias="sidechain_res",
        value_from_simulation_params=True,
        simulation_params_default=30,
        categories=["Simulation parameters"],
        description="Resolution of side chain sampling (degrees).",
    )
    steric_trials: int = Field(
        value_from_simulation_params=True,
        simulation_params_default=250,
        categories=["Simulation parameters"],
        description="Number of tries to get a non-clashing conformation for each step.",
    )
    overlap_factor: float = Field(  # validator between 0 and 1
        value_from_simulation_params=True,
        simulation_params_default=0.65,
        categories=["Simulation parameters"],
        description="Van der Waals overlap factor. Lowering its value will result in higher steric overlap being "
        "accepted.",
    )
    steering: int = Field(
        value_from_simulation_params=True,
        simulation_params_default=0,
        categories=["Simulation parameters"],
        description="Frequency at which steering vector is updated (if using steering).",
    )
    solvent: str = Field(
        value_from_simulation_params=True,
        simulation_params_default="VDGBNP",
        categories=["Simulation parameters"],
    )
    usesrun: bool = Field(
        categories=["General settings"],
        description="Run PELE using SRUN instead of MPI.",
        default=False,
    )
    spawning: str = Field(  # TODO: validator for available options
        value_from_simulation_params="spawning_type",
        simulation_params_default="independent",
        categories=["Simulation parameters"],
        description="Spawning type decides how trajectories are distributed throughout available CPUs depending on the"
        " features of the clusters.",
    )
    iterations: int = Field(
        tests_value=1,
        value_from_simulation_params=True,
        categories=["Simulation parameters"],
        description="Number of iterations to perform.",
    )
    pele_steps: int = Field(
        alias="steps",
        value_from_simulation_params=True,
        simulation_params_default=8,
        tests_value=1,
        categories=["Simulation parameters"],
        description="Number of PELE steps in each iteration.",
    )
    cpus: int = Field(
        tests_value=5,
        value_from_simulation_params=True,
        simulation_params_default=60,
        categories=["Simulation parameters"],
        description="Number of CPUs to use.",
    )
    density: str = Field(  # TODO: validator to make sure user chooses from available options
        value_from_simulation_params=True,
        simulation_params_default="null",
        categories=["Simulation parameters"],
        description="Type of density calculator when running and adaptive simulation. Available options: 'constant' "
        "(default), 'heaviside', 'continuous', 'exitContinuous'.",
    )
    cluster_values: List[float] = Field(
        value_from_simulation_params=True,
        simulation_params_default=[1.75, 2.5, 4, 6],
        categories=["Simulation parameters"],
    )
    cluster_conditions: List[float] = Field(
        value_from_simulation_params=True,
        simulation_params_default=[1, 0.6, 0.4, 0.0],
        categories=["Simulation parameters"],
    )
    equilibration: bool = Field(
        categories=["Simulation parameters"],
        default=False,
        description="Whether to run equilibration or not.",
    )
    clust_type: str = Field(
        value_from_simulation_params=True,
        simulation_params_default="rmsd",
        categories=["Simulation parameters"],
    )
    equil_steps: int = Field(
        alias="equilibration_steps",
        categories=["Simulation parameters"],
        description="Number of equilibration steps to perform.",
    )
    adaptive_restart: bool = Field(
        description="Restart adaptive simulation from the last epoch.",
        categories=["General settings"],
        default=False,
    )
    report_name: str = Field(
        alias="report",
        default="report",
        categories=["Simulation parameters"],
        candidate_for_deprecation=True,
        description="Assign custom names to report file.",
    )
    traj_name: str = Field(
        alias="traj",
        default="trajectory.pdb",
        categories=["Simulation parameters"],
        description="Type of trajectory to save, choose 'trajectory.xtc' or 'trajectory.pdb' (default).",
    )
    epsilon: float = Field(
        value_from_simulation_params=True,
        simulation_params_default=0,
        categories=["Simulation parameters"],
        description="The fraction of the processors that will be assigned according to the selected metric when using "
        "epsilon spawning.",
    )
    out_in: bool = Field(categories=["Out in"], description="Launch OutIn package.")
    bias_column: int = Field(
        value_from_simulation_params=True,
        simulation_params_default=5,
        categories=["Simulation parameters"],
        description="The number of the report column towards which the bias should be applied.",
    )
    gridres: int = Field(
        default=10,
        categories=["Ligand preparation"],
        description="Resolution of rotamer sampling when parametrizing hetero molecules.",
    )
    use_peleffy: bool = Field(
        categories=["Ligand preparation"],
        can_be_falsy=True,
        default=False,
        description="Use peleffy which supports OpenFF forcefield as well as OPLS2005.",
    )
    core: Union[int, List[str]] = Field(
        categories=["Ligand preparation"],
        description="List of PDB atom names that will be included as part of the rigid core. In case it is not "
        "specified, the algorithm will pick up a set of non-rotatable atoms centered in the molecular structure.",
    )
    mtor: int = Field(
        alias="maxtorsion",
        default=4,
        categories=["Ligand preparation"],
        description="Maximum number of rotamers for each flexible ligand branch.",
    )
    n: int = Field(
        default=10000,
        categories=["Ligand preparation"],
        description="Maximum number of flexible side chains in the ligand.",
    )
    solvent_template: str = Field(
        categories=["Ligand preparation"],
        description="Path to file with solvent template.",
    )
    mae_lig: str = Field(
        categories=["Ligand preparation"],
        description="Maestro file containing quantum charges for the ligand preparation.",
    )
    no_ppp: bool = Field(
        alias="skip_preprocess",
        value_from_simulation_params=True,
        simulation_params_default=False,
        categories=["Protein preparation"],
        description="Skip protein preprocessing.",
    )
    skip_prep: bool = Field(
        value_from_simulation_params=True,
        simulation_params_default=False,
        categories=["Ligand preparation"],
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
    nonstandard: List[str] = Field(
        value_from_simulation_params=True,
        simulation_params_default=[],
        categories=["Protein preparation"],
        description="List of non-standard residue names to be omitted during protein preprocessing.",
    )
    prepwizard: bool = Field(
        default=False,
        description="Whether to run Protein Preparation Wizard before launching the simulation.",
        categories=["Protein preparation"],
    )
    box_radius: float = Field(
        value_from_simulation_params=True,
        categories=["Box settings"],
        description="Radius of the simulation box.",
    )
    box_center: Union[List[float], str] = Field(
        categories=["Box settings"],
        description="Center of the simulation box. It can be set using x, y, z coordinates or by specifying a center "
        "atom following the format: 'chain ID:residue number:atom name'.",
    )
    native: str = Field(
        alias="rmsd_pdb",
        categories=["Metrics"],
        description="Path to the PDB file with a reference structure. Mind that both proteins need to be prealigned as "
        "well as share the same chain ID and PDB atom names of the ligand.",
    )
    atom_dist: List[str] = Field(
        default_factory=list,
        categories=["Metrics"],
        description="List of strings representing atoms, between which the distance will be calculated throughout the "
        "simulation.",
    )
    debug: bool = Field(
        default=False,
        categories=["General settings"],
        description="Launch debug mode which generates all input files but does not start the simulation.",
    )
    folder: str = Field(
        alias="working_folder",
        description="Custom name of the top level directory where the simulation will be saved.",
    )
    output: str = Field(
        default="output",
        categories=["General settings"],
        description="Name of the raw output folder.",
    )
    randomize: bool = Field(
        value_from_simulation_params=True,
        simulation_params_default=False,
        categories=["Global Exploration", "Simulation parameters"],
        description="Whether to randomize the initial positions of the ligand around the protein (or center of "
        "interface).",
    )
    # TODO: Check out all exit keywords and what they do.
    exit: bool = Field(categories=["In Out"])
    exit_value: float = Field(
        value_from_simulation_params=True,
        simulation_params_default=0.9,
        categories=["In Out"],
    )
    exit_condition: str = Field(
        categories=["In Out"],
        value_from_simulation_params=True,
        simulation_params_default=">",
    )
    exit_trajnum: int = Field(
        categories=["In Out"],
        value_from_simulation_params=True,
        simulation_params_default=4,
    )
    poses: int = Field(
        description="Number of ligand poses to generate during randomization.",
        value_from_simulation_params=True,
    )
    full: bool = Field(
        alias="global",
        categories=["Global Exploration"],
        description="Launch global exploration package.",
    )
    clust: int = Field(
        alias="exit_clust", tests_value=2, candidate_for_deprecation=True
    )
    restart: bool = Field(
        categories=["General settings"],
        description="Restart the simulation from scratch (not from the last epoch), without "
        "overwriting configuration files.",
    )
    rescoring: bool = Field(
        categories=["Rescoring"], description="Launch rescoring package."
    )
    in_out: bool = Field(categories=["In Out"])
    in_out_soft: bool = Field(categories=["In Out"])
    waters: Union[str, List[str]] = Field(
        value_from_simulation_params=True,
        simulation_params_default=[],
        categories=["Water"],
        description="List of water molecules to perturb which can be defined using residue identifiers (e.g. 'A:123') "
        "or using keyword 'all' to select all of them.",
    )
    water_center: Union[List[float], str] = Field(
        categories=["Water"],
        description="Center of the water box, can be defined as x, y, z coordinates or by selecting an atom as the "
        "center (format 'chain ID:residue number:atom name'.",
    )
    water_temp: float = Field(
        value_from_simulation_params=True,
        simulation_params_default=5000,
        categories=["Water"],
        description="Temperature of water perturbation step.",
    )
    water_overlap: float = Field(
        value_from_simulation_params=True,
        simulation_params_default=0.78,
        categories=["Water"],
        description="Overlap factor for water perturbation step.",
    )
    water_constr: float = Field(
        value_from_simulation_params=True,
        simulation_params_default=0,
        categories=["Water"],
        description="Spring constant of the center of mass constraint applied to the water molecules after "
        "perturbation.",
    )
    water_trials: int = Field(
        value_from_simulation_params=True,
        simulation_params_default=10000,
        categories=["Water"],
        description="Number of trials while performing water perturbation.",
    )
    water_radius: float = Field(
        default=6.0, categories=["Water"], description="Radius of the water box."
    )
    induced_fit_exhaustive: bool = Field(
        categories=["Induced fit"], description="Launch Induced Fit Exhaustive package."
    )
    induced_fit_fast: bool = Field(
        categories=["Induced fit"], description="Launch Induced Fit Fast package."
    )
    ca_constr: float = Field(
        value_from_simulation_params=True,
        simulation_params_default=0.5,
        categories=["Constraints"],
        description="Spring constant for constraining alpha carbons in the backbone.",
    )
    ca_interval: int = Field(
        value_from_simulation_params=True,
        simulation_params_default=10,
        categories=["Constraints"],
        description="Interval at which alpha carbons in the backbone will be constrained (every n atoms).",
    )
    seed: int = Field(
        default_factory=generate_random_seed,
        tests_value=12345,
        categories=["General settings"],
    )
    log: bool = Field(
        default=False,
        description="Retrieve PELE log files during simulation.",
        categories=["General settings"],
    )
    pele_exec: str = Field()
    pele_data: str = Field()
    pele_documents: str = Field(default=os.path.join(constants.PELE, "Documents"))
    pca: str = Field(
        description="Path to PCA file", categories=["PCA"]
    )  # TODO: Difference between this and pca_traj?
    anm_direction: str = Field(
        value_from_simulation_params=True,
        simulation_params_default="random",
        categories=["ANM"],
    )
    anm_mix_modes: str = Field(
        value_from_simulation_params=True,
        simulation_params_default="mixMainModeWithOthersModes",
        categories=["ANM"],
    )
    anm_picking_mode: str = Field(
        value_from_simulation_params=True,
        simulation_params_default="RANDOM_MODE",
        categories=["ANM"],
    )
    anm_displacement: float = Field(
        value_from_simulation_params=True,
        simulation_params_default=0.75,
        categories=["ANM"],
    )
    anm_modes_change: int = Field(
        value_from_simulation_params=True,
        simulation_params_default=4,
        categories=["ANM"],
    )
    anm_num_of_modes: int = Field(
        value_from_simulation_params=True,
        simulation_params_default=6,
        categories=["ANM"],
    )
    anm_relaxation_constr: float = Field(
        value_from_simulation_params=True,
        simulation_params_default=0.5,
        categories=["ANM"],
    )
    remove_constraints: bool = Field(
        can_be_falsy=True,
        value_from_simulation_params=True,
        simulation_params_default=False,
        categories=["ANM"],
    )
    pca_traj: Union[List[str], str] = Field(
        categories=["PCA"],
        description="List of protein PDB files to be used for " "PCA analysis.",
    )
    perturbation: bool = Field(categories=["Advanced"])
    binding_energy: str = Field(
        value_from_simulation_params=True,
        simulation_params_default=constants.BE,
        categories=["Advanced"],
    )
    analyse: bool = Field(
        default=True,
        description="Whether to run analysis or not.",
        categories=["Analysis"],
    )
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
        value_from_simulation_params=True,
        description="Set to 'min' or 'max' to either minimise ot maximise epsilon.",
    )
    external_constraints: List[str] = Field(
        default=[],
        categories=["Constraints"],
        description="List of user-defined positional or atom-atom constraints.",
    )
    only_analysis: bool = Field(
        default=False,
        categories=["Analysis"],
        description="Analyse simulation output without relaunching PELE.",
    )
    analysis_nclust: int = Field(
        default=10,
        categories=["Analysis"],
        description="Number of clusters to extract when using the Gaussian Mixture Model.",
    )
    te_column: int = Field(
        default=4,
        categories=["Analysis"],
        description="Number of report column containing total energy.",
    )
    be_column: int = Field(
        default=5,
        categories=["Analysis"],
        description="Number of report column containing binding energy.",
    )
    limit_column: int = Field(
        default=6,
        categories=["Analysis"],
        description="Number of report column where your the metrics start.",
    )
    com: float = Field(
        alias="COMligandConstraint",
        value_from_simulation_params="COMligandConstraint",
        simulation_params_default=0,
        categories=["FragPELE"],
    )
    pele_license: str = Field(
        default=os.path.join(constants.PELE, "licenses"),
        categories=["General settings"],
    )
    license: str = Field(value_from="pele_license", categories=["General settings"])
    schrodinger: str = Field(
        categories=["General settings"],
        description="Path to Schrodinger, e.g. '/home/software/schrodinger2021-3/'",
    )
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
    permissive_metal_constr: bool = Field(categories=["Constraints"])
    constrain_all_metals: bool = Field(default=False, categories=["Constraints"])
    no_metal_constraints: bool = Field(default=False, categories=["Constraints"])
    frag_run: bool = Field(default=True, categories=["FragPELE"])
    frag_core: str = Field(categories=["FragPELE"])
    frag_input: str = Field(default=None, categories=["FragPELE"])
    frag_ligands: str = Field(default=None, categories=["FragPELE"])
    growing_steps: int = Field(default=6, categories=["FragPELE"])
    frag_steps: int = Field(alias="steps_in_gs", default=3, categories=["FragPELE"])
    frag_eq_steps: int = Field(
        alias="sampling_steps", default=20, categories=["FragPELE"]
    )
    protocol: str = Field(categories=["FragPELE"])
    frag_ai: bool = Field(categories=["FragPELE"])
    frag_ai_iterations: int = Field(categories=["FragPELE"])
    chain_core: str = Field(default="L", categories=["FragPELE"])
    frag_restart: str = Field(categories=["FragPELE"])
    frag_criteria: str = Field(default="Binding Energy", categories=["FragPELE"])
    frag_output_folder: str = Field(default="growing_steps", categories=["FragPELE"])
    frag_cluster_folder: str = Field(default="clustering_PDBs", categories=["FragPELE"])
    frag_library: str = Field(default=None, categories=["FragPELE"])
    frag_core_atom: Union[str, None] = Field(categories=["FragPELE"], default=None)
    analysis_to_point: Optional[List[float]] = Field(categories=["FragPELE"])
    fragment_atom: str = Field(default=None, categories=["FragPELE"])
    frag_restart_libraries: bool = Field(default=None, categories=["FragPELE"])
    proximityDetection: bool = Field()
    ppi: bool = Field(categories=["PPI"])
    center_of_interface: str = Field(
        categories=["PPI"],
        description="Atom string defining the center of interface in PPI simulation. Should "
        "follow the format: 'chain ID:residue number:atom name'.",
    )
    protein: str = Field(categories=["PPI"])
    ligand_pdb: str = Field(
        categories=["PPI"], description="Path to the ligand PDB file in PPI simulation."
    )
    skip_refinement: bool = Field(default=False, categories=["PPI", "Site finder"])
    n_waters: int = Field(
        value_from_simulation_params=True,
        simulation_params_default=0,
        categories=["Water"],
    )
    site_finder: bool = Field(
        categories=["Site finder"], description="Launch site finder package."
    )
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
    charge_parametrization_method: str = Field(
        categories=["Ligand preparation"],
        description="Ligand charge parametrization method when using peleffy. Choose one of: 'gasteiger', 'am1bcc', "
        "'OPLS2005'.",
    )
    exclude_terminal_rotamers: bool = Field(
        default=True, categories=["Ligand preparation"]
    )
    singularity_exec: str = Field(categories=["General settings"])
    terminal_constr: float = Field(
        value_from_simulation_params="terminal_constr", simulation_params_default=5.0
    )
    covalent_residue: str = Field(categories=["Covalent docking"])
    nonbonding_radius: float = Field(
        categories=["Covalent docking"],
        value_from_simulation_params=True,
        simulation_params_default=20.0,
    )
    perturbation_trials: int = Field(
        categories=["Covalent docking"],
        value_from_simulation_params=True,
        simulation_params_default=10,
        description="",
    )
    refinement_angle: float = Field(
        categories=["Covalent docking"],
        value_from_simulation_params=True,
        simulation_params_default=10.0,
        description="Extend of side chain perturbation during covalent docking refinement.",
    )
    covalent_docking_refinement: bool = Field(
        categories=["Covalent docking"],
        description="Launch covalent docking refinement.",
    )
    ligand_conformations: str = Field(categories=["Ligand conformations"])
    conformation_freq: int = Field(
        categories=["Ligand conformations"],
        value_from_simulation_params=True,
        default=4,
    )
    overlap_factor_conformation: float = Field(
        categories=["Ligand conformations"],
        description="Van der Waals overlap factor in conformation perturbation.",
    )
    inter_step_logger: bool = Field(
        categories=["General settings"],
        description="Activate interstep logger.",
        default=False,
    )
    minimum_steps: bool = Field(
        categories=["General settings"],
        description="Force explorers that completed their steps earlier to keep on running PELE steps until all "
        "explorers have finished.",
    )
    site_finder_global: bool = Field(
        categories=["Site finder"], description="Launch site finder global exploration."
    )
    site_finder_local: bool = Field(
        categories=["Site finder"], description="Launch site finder local exploration."
    )
    kde: bool = Field(
        categories=["Analysis"],
        description="Generate kernel density estimate (KDE) plot when analysing the simulation.",
        default=False,
    )
    kde_structs: int = Field(
        categories=["Analysis"],
        default=1000,
        description="Number of structures to include on the KDE plot.",
    )
    plot_filtering_threshold: float = Field(
        categories=["Analysis"],
        description="Percentage of output structures to filter out before creating plots.",
    )
    clustering_filtering_threshold: float = Field(
        categories=["Analysis"],
        default=0.25,
        can_be_falsy=True,
        description="Percentage of output structures to filter our before clustering",
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
        description="Method for selecting representative structures for each cluster, "
        "you can choose one of: 'total_25_percentile', 'total_5_percentile', "
        "'total_mean', 'interaction_25_percentile', 'interaction_5_percentile' or 'interaction_mean'.",
    )
    bandwidth: float = Field(
        categories=["Analysis"],
        value_from_simulation_params=True,
        simulation_params_default=2.5,
        description="Bandwidth for the mean shift and HDBSCAN clustering (also called epsilon).",
    )
    max_top_clusters: int = Field(
        categories=["Analysis"],
        can_be_falsy=True,
        default=8,
        description="Maximum number of clusters to retrieve during simulation analysis.",
    )
    min_population: float = Field(
        categories=["Analysis"],
        value_from_simulation_params=True,
        simulation_params_default=0.01,
        can_be_falsy=True,
        description="The minimum amount of structures in each cluster, takes a value between 0 and 1, "
        "where 0.01 refers to 1%.",
    )
    max_top_poses: int = Field(
        categories=["Analysis"],
        can_be_falsy=True,
        default=100,
        description="The maximum number of top poses (structures with lowest binding energy) to be "
        "retrieved during simulation analysis.",
    )
    saturated_mutagenesis: bool = Field(
        categories=["Saturated mutagenesis"],
        description="Launch saturated mutagenesis package.",
    )
    cpus_per_mutation: int = Field(
        categories=["Saturated mutagenesis"],
        tests_value=2,
        description="Number of CPUs to use for each mutation.",
    )
    constraint_level: int = Field(
        categories=["Constraints"],
        description="Select constraint level with predefined parameters parameter. Accepts value from 0 (no "
        "constraints) to 3.",
        value_from_simulation_params=True,
    )
    mae: bool = Field(categories=["Frag"], description="Save FragPELE output as MAE files.")

    @validator("*", pre=True, always=True, allow_reuse=True)
    def set_tests_values(cls, v, values, field):
        """
        Adjusts values of parameters when running tests.
        """
        if values.get("test"):
            test_value = field.field_info.extra.get("tests_value")
            if test_value is not None:
                return test_value
        return v

    @validator("*", pre=True, always=True, allow_reuse=True)
    def set_value_from(cls, v, values, field):
        """
        Gets value from an a duplicated parameter, e.g. steps and pele_steps.
        """
        value_from = field.field_info.extra.get("value_from")
        if value_from is not None:
            return values.get(value_from)
        return v

    @validator("system", "mae_lig", allow_reuse=True)
    def construct_path(cls, v):
        """
        Gets absolute path to the file.
        """
        return str(os.path.abspath(v)) if v else None

    @validator("residue", allow_reuse=True)
    def validate_residue_name(cls, v):
        """
        Checks the residue name for unsupported 'UNK' value.
        """
        if v == "UNK":
            raise custom_errors.LigandNameNotSupported(
                "'UNK' ligand name is not supported, please rename it, e.g. 'LIG'."
            )
        return v

    @validator("usesrun", always=True, allow_reuse=True)
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
        allow_reuse=True,
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

    @validator("covalent_residue", allow_reuse=True)
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

    @validator("frag_core_atom", allow_reuse=True)
    def validate_frag_core_atom(cls, v):
        """
        Checks if the string fits a regex pattern indicating the core atom,e.g. C3-H2.
        """
        pattern = r"([A-Z]{1,2}\d{1,2}\-H[A-Z]?\d{1,2})"

        if v:
            if not re.match(pattern, v):
                raise custom_errors.WrongAtomStringFormat(
                    f"Atom string set in {v} does not seem to have the right format. It should follow C atom name: H "
                    "atom name format."
                )
        return v

    @validator("sidechain_resolution", "gridres", allow_reuse=True)
    def check_divisibility(cls, v):
        if v and not 360 % v == 0:
            raise ValueError(
                "The value should be easily multiplied to obtain 360, e.g. 30, 45 or 10 would be valid."
            )
        return v

    @validator("constraint_level", allow_reuse=True)
    def parse_constraint_level(cls, v):
        if v:
            if v not in (0, 1, 2, 3):
                raise ValueError(
                    f"Invalid constraint level {v}, should be 0, 1, 2 or 3."
                )
        return v

    @validator("terminal_constr", "ca_constr", allow_reuse=True)
    def assert_positive_number(cls, v):
        if v and v < 0:
            raise ValueError(f"Value {v} should be positive.")
        return v

    @validator(
        "epsilon", "min_population", "plot_filtering_threshold", allow_reuse=True
    )
    def assert_value_0_1(cls, v):
        if v and not 0 <= v <= 1:
            raise ValueError(f"Value {v} should be between 0 and 1.")
        return v

    @validator("cpus_per_mutation", allow_reuse=True)
    def check_cpus(cls, v, values):
        if v and v >= values.get("cpus"):
            raise ValueError(
                f"The number of cpus_per_mutation cannot be higher than total cpus - 1."
            )
        return v

    @validator(
        "cluster_representatives_criterion",
        "top_clusters_criterion",
        "clustering_method",
        allow_reuse=True,
    )
    def check_cluster_representatives_criterion(cls, v, field):
        """
        Checks if analysis flags are set to one of the available values.
        """
        if v:
            available_values = eval(f"constants.{field.name}_values")

            if isinstance(available_values, dict):
                available_values = available_values.keys()

            if v.lower() not in available_values:
                raise ValueError(
                    f"Value of {field.name} should be one of: {available_values}"
                )
        return v

    @validator("charge_parametrization_method")
    def set_charge_parametrization_method(cls, v, values):
        # TODO
        pass

    @validator("external_constraints")
    def parse_external_constraints(cls, v, values):
        if v:
            return helpers.retrieve_constraints_for_pele(v, values["system"])
        else:
            return []
