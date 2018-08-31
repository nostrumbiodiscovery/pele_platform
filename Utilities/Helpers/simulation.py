import os
from PELE_Platform.Utilities.Helpers import helpers, template_builder
import PELE_Platform.AdaptivePELE.adaptiveSampling as ad
import PELE_Platform.constants as cs
import PELE_Platform.Utilities.Helpers.center_of_mass as cm 

class SimulationBuilder(template_builder.TemplateBuilder):

    def __init__(self, adaptive, pele, env):
        self.adaptive_file = adaptive
        self.pele_file = pele
        self.topology = env.topology
        print(self.adaptive_file, self.pele_file, env.__dict__)
        self.fill_pele_template(env)
        self.fill_adaptive_template(env)

    def fill_pele_template(self, env):
        self.pele_keywords = { "NATIVE": env.native, "FORCEFIELD": env.forcefield, "CHAIN": env.chain, 
                        "CONSTRAINTS": "\n".join(env.constraints), "CPUS":env.cpus,
                        "LICENSES": cs.LICENSE, "BOX_RADIUS": env.box_radius, "BOX_CENTER": env.box_center, "HBOND1": env.hbond_donor, 
                        "HBOND2": env.hbond_acceptor, "SASA_min": env.sasa_min, "SASA_max": env.sasa_max }
        super(SimulationBuilder, self).__init__(self.pele_file, self.pele_keywords)

    def fill_adaptive_template(self, env):
        self.adaptive_keywords = { "RESTART": cs.RESTART, "OUTPUT": env.adap_ex_output, "INPUT":env.adap_ex_input,
                "CPUS":env.cpus, "PELE_CFILE": self.pele_file, "LIG_RES": env.residue, "SEED": env.seed, "EQ_STEPS": env.equil_steps}
        super(SimulationBuilder, self).__init__(self.adaptive_file, self.adaptive_keywords)



    def run(self, hook=False):
        with helpers.cd(os.path.dirname(self.file)):
            if hook:
		ad.main(self.file, clusteringHook=self.interactive_clustering)
            else:
		ad.main(self.file)

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


