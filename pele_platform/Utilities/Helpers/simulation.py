from dataclasses import dataclass
from copy import deepcopy
import os

import AdaptivePELE.adaptiveSampling as ad
from pele_platform.Utilities.Helpers import helpers, template_builder
import pele_platform.Utilities.Parameters.parameters as pv
import pele_platform.Utilities.Helpers.water as wt
from pele_platform.Utilities.Parameters.parameters import Parameters
from pele_platform.constants import constants
from pele_platform.Adaptive.interaction_restrictions import InteractionRestrictionsBuilder
from pele_platform.context import context


@dataclass
class SimulationBuilder(template_builder.TemplateBuilder):
    adaptive_file: str
    pele_file: str
    topology: str

    def generate_inputs(self, water_obj):
        """
        It generates the input files for Adaptive PELE, according to the
        parameters previously generated.

        Parameters
        ----------
        water_obj : a WaterIncluder object
            The parameters for aquaPELE, if applicable
        """
        # Convert raw parameters to JSON strings
        formatted_params = self.format_parameters()

        # Fill two times because we have flags inside flags
        self.fill_pele_template(formatted_params, water_obj)
        self.fill_pele_template(formatted_params, water_obj)
        self.fill_adaptive_template(formatted_params)
        self.fill_adaptive_template(formatted_params)

    @staticmethod
    def format_parameters():
        """
        Injects user-defined parameters into strings compatible with PELE configuration files.

        Returns
        --------
        Deepcopy of parameters with correct JSON string formatting.
        """
        parameters = context.parameters_copy

        if parameters.covalent_residue:
            parameters.covalent_sasa = constants.SASA_COVALENT.format(parameters.covalent_residue)
            parameters.max_trials_for_one = parameters.perturbation_trials * 2
            parameters.sidechain_perturbation = constants.SIDECHAIN_PERTURBATION

            if parameters.refinement_angle and parameters.covalent_docking_refinement:
                parameters.refinement_angle = constants.refinement_angle.format(parameters.refinement_angle)
            else:
                parameters.refinement_angle = ""

        else:
            parameters.refinement_angle = ""
            parameters.covalent_sasa = ""
            parameters.sidechain_perturbation = ""

        if not parameters.interaction_restrictions:
            parameters.met_interaction_restrictions = ""
            parameters.interaction_restrictions = ""
        else:
            ir_parser = InteractionRestrictionsBuilder()
            ir_parser.parse_interaction_restrictions(parameters.system, parameters.interaction_restrictions)
            parameters.parameters = ir_parser.fill_template(parameters.parameters)
            parameters.met_interaction_restrictions = ir_parser.metrics_to_json()
            parameters.interaction_restrictions = ir_parser.conditions_to_json()

        if parameters.pca:
            parameters.pca = constants.PCA.format(parameters.pca)
        else:
            parameters.pca = ""

        if parameters.singularity_exec:
            parameters.mpi_params = parameters.singularity_exec
            if parameters.frag_pele:
                parameters.pele_exec = parameters.mpi_params + " Pele_mpi"
            else:
                parameters.pele_exec = "Pele_mpi"

        parameters.water = parameters.water if parameters.water else ""
        parameters.allow_empty_selectors = '"allowEmptyWaterSelectors": true,' if parameters.allow_empty_selectors else ""

        mpi_params_name = "srunParameters" if parameters.usesrun else "mpiParameters"
        parameters.mpi_params = (
            f'"{mpi_params_name}": "{parameters.mpi_params}",' if parameters.mpi_params else ""
        )

        parameters.usesrun = str(parameters.usesrun).lower()
        parameters.minimum_steps = constants.MINIMUMSTEPS if parameters.minimum_steps else ""
        parameters.inter_step_logger = constants.INTERSTEPLOGGER if parameters.inter_step_logger else ""
        parameters.conformation_perturbation = constants.CONFORMATION_PERTURBATION if parameters.ligand_conformations else ""
        parameters.verbose = str(parameters.verbose).lower()
        parameters.equilibration = str(parameters.equilibration).lower()

        parameters.logfile = parameters.logfile if parameters.logfile else '"simulationLogPath" : "$OUTPUT_PATH/logFile.txt",'
        parameters.spawning_condition = f'"condition": "{parameters.spawning_condition}",' if parameters.spawning_condition else ""
        parameters.proximityDetection = str(parameters.proximityDetection).lower()

        return parameters

    def fill_pele_template(
            self, env: Parameters, water_obj: wt.WaterIncluder
    ) -> None:
        # Translate OpenFF force field names into the format expected by PELE
        if "openff" in env.forcefield.lower():
            forcefield = "OpenFF-OPLS2005"
        else:
            forcefield = "OPLS2005"

        # Fill in PELE template
        self.pele_keywords = {
            "PERTURBATION": env.perturbation,
            "SIDECHAIN_PERTURBATION": env.sidechain_perturbation,
            "SELECTION_TO_PERTURB": env.selection_to_perturb,
            "BE": env.binding_energy,
            "SASA": env.sasa,
            "LOGFILE": env.logfile,
            "NATIVE": env.native,
            "FORCEFIELD": forcefield,
            "CHAIN": env.chain,
            "CONSTRAINTS": "\n".join(env.constraints),
            "CPUS": env.cpus,
            "LICENSES": env.license,
            "BOX_RADIUS": env.box_radius,
            "BOX_CENTER": env.box_center,
            "SASA_min": env.sasa_min,
            "SASA_max": env.sasa_max,
            "WATER_RADIUS": water_obj.water_radius,
            "WATER_CENTER": water_obj.water_center,
            "WATER": water_obj.water_line,
            "WATER_ENERGY": water_obj.water_energy,
            "METRICS": env.metrics,
            "REPORT_NAME": env.report_name,
            "TRAJ_NAME": env.traj_name,
            "SOLVENT": env.solvent.upper(),
            "PARAMETERS": env.parameters,
            "SIDECHAIN_RESOLUTION": env.sidechain_resolution,
            "OVERLAP": env.overlap_factor,
            "STERIC_TRIALS": env.steric_trials,
            "TEMPERATURE": env.temperature,
            "MIN_FREQ": env.min_freq,
            "SIDECHAIN_FREQ": env.sidechain_freq,
            "WATER_FREQ": env.water_freq,
            "COMFPERT_FREQ": env.conformation_freq,
            "ANM_FREQ": env.anm_freq,
            "BOX": env.box,
            "PROXIMITY": env.proximityDetection,
            "VERBOSE": env.verbose,
            "ANM_DISPLACEMENT": env.anm_displacement,
            "ANM_MODES_CHANGE": env.anm_modes_change,
            "ANM_DIRECTION": env.anm_direction,
            "ANM_MIX_MODES": env.anm_mix_modes,
            "ANM_PICKING_MODE": env.anm_picking_mode,
            "ANM_NUM_OF_MODES": env.anm_num_of_modes,
            "ANM_RELAXATION_CONST": env.anm_relaxation_constr,
            "PCA": env.pca,
            "COMPLEXES": env.complexes,
            "PELE_STEPS": env.frag_pele_steps,
            "OUTPUT_PATH": env.output_path,
            "COM": env.com,
            "STEERING": env.steering,
            "MET_INTERACTION_RESTRICTIONS": env.met_interaction_restrictions,
            "INTERACTION_RESTRICTIONS": env.interaction_restrictions,
            "INTER_STEP_LOGGER": env.inter_step_logger,
            "COVALENT_RESIDUE": env.covalent_residue,
            "REFINEMENT_ANGLE": env.refinement_angle,
            "LOCAL_NONBONDING_ENERGY": env.local_nonbonding_energy,
            "TRIALS": env.perturbation_trials,
            "COVALENT_SASA": env.covalent_sasa,
            "MAXTRIALSFORONE": env.max_trials_for_one,
            "MINIMUM_STEPS": env.minimum_steps,
            "CONFORMATION_PERTURBATION": env.conformation_perturbation,
            "OVERLAP_CONFORMATION": env.overlap_factor_conformation,
        }

        super(SimulationBuilder, self).__init__(self.pele_file, self.pele_keywords)

    def fill_adaptive_template(self, parameters: pv.Parameters) -> None:
        # Fill in adaptive template
        self.adaptive_keywords = {
            "RESTART": parameters.adaptive_restart,
            "OUTPUT": parameters.output,
            "INPUT": parameters.adap_ex_input,
            "CPUS": parameters.cpus,
            "PELE_CFILE": os.path.basename(self.pele_file),
            "LIG_RES": parameters.residue,
            "SEED": parameters.seed,
            "EQ_STEPS": parameters.equil_steps,
            "EQUILIBRATION": parameters.equilibration,
            "EQUILIBRATION_MODE": parameters.equilibration_mode,
            "EPSILON": parameters.epsilon,
            "BIAS_COLUMN": parameters.bias_column,
            "ITERATIONS": parameters.iterations,
            "PELE_STEPS": parameters.pele_steps,
            "REPORT_NAME": parameters.report_name,
            "SPAWNING_TYPE": parameters.spawning,
            "DENSITY": parameters.density,
            "SIMULATION_TYPE": parameters.simulation_type,
            "CLUSTER_VALUES": parameters.cluster_values,
            "CLUSTER_CONDITION": parameters.cluster_conditions,
            "UNBINDING": parameters.unbinding_block,
            "USESRUN": parameters.usesrun,
            "LIGAND": parameters.ligand,
            "PELE_BIN": parameters.pele_exec,
            "PELE_DATA": parameters.pele_data,
            "PELE_DOCUMENTS": parameters.pele_documents,
            "CONDITION": parameters.spawning_condition,
            "MPI_PARAMS": parameters.mpi_params,
            "CLUST_TYPE": parameters.clust_type,
        }
        super(SimulationBuilder, self).__init__(
            self.adaptive_file, self.adaptive_keywords
        )

    def run(self, hook=False) -> None:
        # Launch montecarlo simulation
        with helpers.cd(os.path.dirname(self.adaptive_file)):
            ad.main(self.adaptive_file)
