import os
import AdaptivePELE.adaptiveSampling as ad
from pele_platform.Utilities.Helpers import helpers, template_builder
import pele_platform.constants.constants as cs
import pele_platform.Utilities.Helpers.center_of_mass as cm 

class SimulationBuilder(template_builder.TemplateBuilder):

    def __init__(self, adaptive, pele, env):
        self.adaptive_file = adaptive
        self.pele_file = pele
        self.topology = env.topology
        self.fill_pele_template(env)
        self.fill_pele_template(env)
        self.fill_adaptive_template(env)
        self.fill_adaptive_template(env)

    def fill_pele_template(self, env):
        self.pele_keywords = { "PERTURBATION": env.perturbation, "SELECTION_TO_PERTURB": env.selection_to_perturb,
                        "BE": env.binding_energy, "SASA": env.sasa,
                        "LOGFILE": env.logfile, "NATIVE": env.native, "FORCEFIELD": env.forcefield, "CHAIN": env.chain, 
                        "CONSTRAINTS": "\n".join(env.constraints), "CPUS":env.cpus,
                        "LICENSES": cs.LICENSE, "BOX_RADIUS": env.box_radius, "BOX_CENTER": env.box_center, "HBOND1": env.hbond_donor, 
                        "HBOND2": env.hbond_acceptor, "SASA_min": env.sasa_min, "SASA_max": env.sasa_max,
                        "WATER_RADIUS": env.water_radius, "WATER_CENTER": env.water_center, "WATER": env.water,
                        "WATER_ENERGY": env.water_energy, "METRICS": env.metrics, "REPORT_NAME": env.report_name, "TRAJ_NAME": env.traj_name,
                        "SOLVENT": env.solvent, "PARAMETERS": env.parameters, "SIDECHAIN_RESOLUTION": env.sidechain_resolution,
                        "OVERLAP": env.overlap_factor, "STERIC_TRIALS": env.steric_trials, "TEMPERATURE": env.temperature, 
                        "MIN_FREQ": env.min_freq, "SIDECHAIN_FREQ": env.sidechain_freq, "ANM_FREQ": env.anm_freq, "BOX" : env.box, "PROXIMITY": env.proximityDetection, "WATER_FREQ": env.water_freq, "VERBOSE": env.verbose, "ANM_DISPLACEMENT": env.anm_displacement, "ANM_MODES_CHANGE": env.anm_modes_change, "ANM_DIRECTION": env.anm_direction, "ANM_MIX_MODES": env.anm_mix_modes, "ANM_PICKING_MODE": env.anm_picking_mode,
                        "ANM_NUM_OF_MODES": env.anm_num_of_modes, "ANM_RELAXATION_CONST": env.anm_relaxation_constr,
                        "PCA": env.pca}


        super(SimulationBuilder, self).__init__(self.pele_file, self.pele_keywords)

    def fill_adaptive_template(self, env):
        self.adaptive_keywords = { "RESTART": env.adaptive_restart, "OUTPUT": env.output, "INPUT":env.adap_ex_input,
                "CPUS":env.cpus, "PELE_CFILE": os.path.basename(self.pele_file), "LIG_RES": env.residue, "SEED": env.seed, "EQ_STEPS": env.equil_steps,
                "EQUILIBRATION":env.equilibration, "EPSILON": env.epsilon, "BIAS_COLUMN": env.bias_column, "ITERATIONS": env.iterations, 
                "PELE_STEPS": env.pele_steps, "REPORT_NAME": env.report_name, "SPAWNING_TYPE": env.spawning, "DENSITY": env.density,
                "SIMULATION_TYPE": env.simulation_type, "CLUSTER_VALUES": env.cluster_values, "CLUSTER_CONDITION": env.cluster_conditions,
                "UNBINDING": env.unbinding_block, "USESRUN": env.usesrun, "LIGAND": env.ligand,
                "PELE_BIN": env.pele_exec, "PELE_DATA": env.pele_data, "PELE_DOCUMENTS": env.pele_documents}
        super(SimulationBuilder, self).__init__(self.adaptive_file, self.adaptive_keywords)



    def run(self, hook=False):
        with helpers.cd(os.path.dirname(self.adaptive_file)):
            if hook:
                ad.main(self.adaptive_file, clusteringHook=self.interactive_clustering)
            else:
                ad.main(self.adaptive_file)

    def interactive_clustering(self, cluster_object, paths, simulationRunner, epoch_number):
        initial_rmsd_cluster_values = cluster_object.thresholdCalculator.values
        while len(cluster_object.clusters.clusters) == 1:
            current_values = cluster_object.thresholdCalculator.values
            cluster_object.thresholdCalculator.values = [ value-0.5 if value > 0.5 else value for value in current_values]
            if cluster_object.thresholdCalculator.values == current_values: break
            cluster_object.emptyClustering()
            ad.clusterPreviousEpochs(cluster_object, epoch_number, paths.epochOutputPathTempletized, simulationRunner, self.topology)
            print("Lowering cluster RMSD to: {}".format(cluster_object.thresholdCalculator.values))
        cluster_object.thresholdCalculator.values = initial_rmsd_cluster_values
        print("Interactive clustering ended restoring initial value {}".format(cluster_object.thresholdCalculator.values))


