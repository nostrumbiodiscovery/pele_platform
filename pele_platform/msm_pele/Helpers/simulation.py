import os
import signal
import time
from msm_pele.Helpers import helpers, template_builder
import msm_pele.constants as cs
import AdaptivePELE.adaptiveSampling as ad
try:
    import subprocess32 as subprocess
except ImportError:
    import subprocess

class SimulationBuilder(template_builder.TemplateBuilder):

    def __init__(self, file, topology, keywords, *args, **kwargs):

        self.file = file
        self.topology = topology
        self.keywords = keywords
        self.ad_opt_new = ["false" if opt is False else opt for opt in locals()["args"]]
        self.ad_opt = ["true" if opt is True else opt for opt in self.ad_opt_new]
        self.replace = {keyword : value for keyword, value in zip(self.keywords, self.ad_opt)}
        self.fill_in_control_file()

    def fill_in_control_file(self):
        super(SimulationBuilder, self).__init__(self.file, self.replace)


    def run_pele(self, env, limitTime=False):
        with helpers.cd(os.path.dirname(self.file)):
            toRun = ["mpirun", "-np", str(env.cpus), cs.PELE_BIN, env.pele_temp]
            print(" ".join(toRun))
            startTime = time.time()

            if limitTime:
                try:
                    proc = subprocess.Popen(toRun, shell=False,  universal_newlines=True, preexec_fn=os.setsid)
                    (out, err) = proc.communicate(timeout=limitTime)
                except subprocess.TimeoutExpired:
                    print("killing")
                    os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
            else:
                proc = subprocess.Popen(toRun, shell=False, universal_newlines=True)
                (out, err) = proc.communicate()
                print(out)
                if err:
                    print(err)

            endTime = time.time()
        return "PELE took %.2f sec" % (endTime - startTime)


    def run_adaptive(self, env, hook=False, limitTime=False):
        with helpers.cd(os.path.dirname(self.file)):
            if hook:
                ad.main(self.file)
            else:
                ad.main(self.file, msm_env=env, limitTime=limitTime)

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
