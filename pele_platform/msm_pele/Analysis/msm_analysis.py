import os
import glob
import shutil
import numpy as np
import msm_pele.Helpers.tica as td
import msm_pele.Helpers.helpers as hp
import msm_pele.Analysis.plotMSMAdvancedInfo as pt
import msm_pele.Analysis.MSM_report as rp
from AdaptivePELE.freeEnergies import extractCoords, prepareMSMFolders, estimateDGAdaptive
from AdaptivePELE.freeEnergies import getRepresentativeStructures as getRepr

# DEFAULT VALUES
LAGTIME = 100
NCLUSTER = 200
CLUSTERINSTRIDE = 1
DIMENSIONS = 2
# select the first run, as it uses all trajectories (no
# bootstrap)
REPRESENTATIVES_FILE = "representative_structures/representative_structures_0.dat"
REPRESENTATIVES_STRUCTURES = "representative_structures_pdb_%d"
PMF_FILE = "pmf_xyzg_0.dat"
# Number of the best structures to write
N_BEST = 5


def analyse_results(env, runTica=True, last=False):
    if env.restart in ["all", "adaptive", "pele", "msm"]:
    	run_msm(env, runTica)
    if env.restart in ["all", "adaptive", "pele", "msm", "analyse"]:
        # In case of more than one simulation, i.e. MSM_0, MSM_1, etc
        for i, folder in enumerate(glob.glob(os.path.join(env.adap_l_output, "MSM_*"))):
            analyse_msm(i, env, folder)
        if last:
            try:
                rp.report_MSM(env, os.path.join(env.adap_l_output, "MSM_{}".format(len(glob.glob(os.path.join(env.adap_l_output, "MSM_*")))-1)))
            except IOError:
                pass


def run_msm(env, runTica=True):
    with hp.cd(env.adap_l_output):
        trajs_per_epoch = len(glob.glob(os.path.join("0", "*report*")))
        if runTica:
            td.main(DIMENSIONS, clusters, env.residue, lagtime, trajs_per_epoch, 1000)
            return()
        else:
            extractCoords.main(lig_resname=env.residue, non_Repeat=False, atom_Ids="", nProcessors=env.cpus, parallelize=False, topology=env.topology)
            prepareMSMFolders.main()
            estimateDGAdaptive.main(trajs_per_epoch, env.lagtime, env.msm_clust, lagtimes=env.lagtimes, output=env.results)
            results_file = summerize(env.results)

def analyse_msm(iteration, env, folder):
    with hp.cd(env.adap_l_output):
        try:
            getRepr.main(os.path.join(env.adap_l_output, folder, REPRESENTATIVES_FILE), ".", output=REPRESENTATIVES_STRUCTURES % iteration, topology=env.topology)
        except IndexError:
            pass
        pt.main(4, 1, 5, folder, True, True, True, None, None, env.system_fix, True, False, None, folder, env.residue)
        

def summerize(pele_path):
    results_file = os.path.join(pele_path, "results.txt")
    with open(results_file, 'r') as results:
        lines = hp.preproces_lines(results.readlines())
        for i, line in enumerate(lines):
            try:
                _, dg, stdDg, _, _ = line
                convergence = asses_convergence(dg, stdDg)
            except ValueError:
                pass
            else:
                line.append(convergence)
            finally:
                lines[i] = " ".join(line)
    with open(results_file, 'w') as results:
        results.write("\n".join(lines)+"\n")
    return results_file


def asses_convergence(dg, stdDg):
    """
       Asses whether the MSM analysis
           was good (G), medium (M) or bad (B).
    """
    convergence = "M"
    convergence_rate = round((abs(float(stdDg)*100)/float(dg)))
    if convergence_rate < 5:
        convergence = "G"
    elif convergence_rate > 10:
        convergence = "B"
    return convergence


def copy_best_structures(pmf_file, output_folder, n_best=5):
    """
        Copy the structures of the states with the lowest pmf value as
        best_x_state_y.pdb, where x is the order of lowest pmf (the absolute
        min is 1, the second 2, etc.) and y is the number of the state
    """
    dest_file = os.path.join(output_folder, "best_%d_state_%d.pdb")
    origin_file = os.path.join(output_folder, "cluster_%d.pdb")
    pmf_values = np.loadtxt(pmf_file)
    # sort by the last column, which corresponds to the pmf
    sorted_pmf = np.argsort(pmf_values[:, -1])
    for i, state in enumerate(sorted_pmf[:n_best]):
        shutil.copy(origin_file % state, dest_file % (i+1, state))


if __name__ == "__main__":
    analyse_results("/home/dsoler/STR_PEle/env.adap_l_output", "STR")
