import sys
import os
from AdaptivePELE.freeEnergies import extractCoords, prepareMSMFolders, estimateDGAdaptive
import MSM_PELE.Helpers.helpers as hp

TRAJS_PER_EPOCH = 50
LAGTIME = 100
NCLUSTER = 200
CLUSTERINSTRIDE = 10


def analyse_results(output_pele, ligand_resname, atom_ids=""):
    with hp.cd(output_pele):
        extractCoords.main(lig_resname=ligand_resname, non_Repeat=True, atom_Ids=atom_ids)
        prepareMSMFolders.main()
        estimateDGAdaptive.main(TRAJS_PER_EPOCH, LAGTIME, NCLUSTER, CLUSTERINSTRIDE)
        summerize(output_pele)


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
        results.write("\n".join(lines))

def asses_convergence(dg, stdDg):
    """
       Asses whether the MSM analysis
           was good (G), medium (M) or bad (B).
    """
    convergence = "M"   
    convergence_rate = round((abs(float(stdDg)*100)/float(dg)))
    if(convergence_rate < 5):
        convergence = "G"
    elif(convergence_rate > 10):
        convergence = "B"
    return convergence
        
                    
    





if __name__ == "__main__":
        analyse_results("/home/dsoler/STR_PEle/output_pele", "STR")
