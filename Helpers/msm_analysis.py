import os
import glob
from AdaptivePELE.freeEnergies import extractCoords, prepareMSMFolders, estimateDGAdaptive
from AdaptivePELE.freeEnergies import getRepresentativeStructures as getRepr
import MSM_PELE.Helpers.helpers as hp
import shutil
import numpy as np

# DEFAULT VALUES
LAGTIME = 100
NCLUSTER = 200
CLUSTERINSTRIDE = 1
# select the first run, as it uses all trajectories (no
# bootstrap)
REPRESENTATIVES_FILE = "representative_structures/representative_structures_0.dat"
REPRESENTATIVES_STRUCTURES = "representative_structures_pdb_%d"
PMF_FILE = "pmf_xyzg_0.dat"
# Number of the best structures to write
N_BEST = 5


def analyse_results(output_pele, ligand_resname, cpus, pele_dir, atom_ids=""):
    with hp.cd(output_pele):
        trajs_per_epoch = len(glob.glob(os.path.join("0", "*traj*")))
        extractCoords.main(lig_resname=ligand_resname, non_Repeat=False, atom_Ids=atom_ids, parallelize=False, nProcessors=cpus)
        prepareMSMFolders.main()
        estimateDGAdaptive.main(trajs_per_epoch, LAGTIME, NCLUSTER, CLUSTERINSTRIDE)
        results_file = summerize(output_pele)
        shutil.move(results_file, os.path.join(pele_dir, "results.txt"))
        # In case of more than one simulation, i.e. MSM_0, MSM_1, etc
        MSM_folders = glob.glob(os.path.join(output_pele, "MSM_*"))
        for i, folder in enumerate(MSM_folders):
            getRepr.main(os.path.join(output_pele, folder, REPRESENTATIVES_FILE), ".", output=REPRESENTATIVES_STRUCTURES % i)
            copy_best_structures(os.path.join(output_pele, folder, PMF_FILE), REPRESENTATIVES_STRUCTURES % i, n_best=N_BEST)


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
    analyse_results("/home/dsoler/STR_PEle/output_pele", "STR")
