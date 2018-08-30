import os
import glob
from MSM_PELE.AdaptivePELE.freeEnergies import extractCoords, prepareMSMFolders, estimateDGAdaptive
from MSM_PELE.AdaptivePELE.freeEnergies import getRepresentativeStructures as getRepr
import MSM_PELE.Helpers.tica as td
import MSM_PELE.Helpers.helpers as hp
import shutil
import numpy as np

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


def analyse_results(env, args, runTica=True):
    lagtime = 1 if args.test else env.lagtime
    lagtimes = None if args.test else None
    clusters = 2 if args.test else env.msm_clust
    with hp.cd(env.adap_l_output):
    	trajs_per_epoch = len(glob.glob(os.path.join("0", "*report*")))
        if runTica:
            td.main(DIMENSIONS, clusters, args.residue, lagtime, trajs_per_epoch, 1000)
            return()
        else:
            extractCoords.main(lig_resname=args.residue, non_Repeat=False, atom_Ids="", nProcessors=args.cpus, parallelize=True, topology=env.topology)
            prepareMSMFolders.main()
            estimateDGAdaptive.main(trajs_per_epoch, lagtime, clusters, lagtimes=lagtimes)
            results_file = summerize(env.adap_l_output)
            shutil.move(results_file, os.path.join(env.pele_dir, "results.txt"))
            # In case of more than one simulation, i.e. MSM_0, MSM_1, etc
            MSM_folders = glob.glob(os.path.join(env.adap_l_output, "MSM_*"))
            for i, folder in enumerate(MSM_folders):
		try:
		    getRepr.main(os.path.join(env.adap_l_output, folder, REPRESENTATIVES_FILE), ".", output=REPRESENTATIVES_STRUCTURES % i, topology=env.topology)
		except IndexError: 
		    pass

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
