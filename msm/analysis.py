import sys
import os
from AdaptivePELE.freeEnergies import extractCoords, prepareMSMFolders, estimateDGAdaptive
sys.path.append("/home/dsoler/PelePlop")
import Helpers.pele_env as pele

TRAJS_PER_EPOCH = 50
LAGTIME = 100
NCLUSTER = 200
CLUSTERINSTRIDE = 10


def analyse_results(output_pele, ligand_resname, atom_ids=""):
	with pele.cd(output_pele):
		extractCoords.main(lig_resname=ligand_resname, non_Repeat=True, atom_Ids=atom_ids)
		prepareMSMFolders.main()
		estimateDGAdaptive.main(TRAJS_PER_EPOCH, LAGTIME, NCLUSTER, CLUSTERINSTRIDE)



def summerize(pele_paths, residues):
	gibs_energies  = ["DG LIGAND RANKING", "-------------"]
	for pele_path in pele_paths:
		results_file = os.path.join(pele_path, "output_adaptive_long/results.txt")
	        with open(results_file, 'r') as results:
			gibs_energies = [energy.strip("\n") for energy in results if not energy.startswith("#")]
	results = [ [gib_energy, residue] for gib_energy, residue in zip(gibs_energies, residues)]
	results.sort(key=lambda x: x[0])
        convergences = asses_convergence(results)
	final_report = ["{}: {} {}".format(result[1], result[0], conv) for (result, conv) in zip(results, convergences)]
	return final_report

def asses_convergence(results):
	"""
	   Asses whether the MSM analysis
           was good (G), medium (M) or bad (B).
	"""
	convergence = "M"
        convergences = []
	
	for result in results:
            values, _ = result
            epoch, dg, stdDg, db, StdD = values.split()
            convergence_rate = round((abs(float(stdDg)*100)/float(dg)))
            if(convergence_rate < 5):
	 	convergence = "G"
            elif(convergence_rate > 10):
		convergence = "B"
            convergences.append(convergence)
	return convergences
		
					
	





if __name__ == "__main__":
	output=summerize(["/home/dsoler/STR_PEle",], ["STR",])
        output.insert(0,"#Resiude Epoch DG StdDG Db StdDb\n#==============================\n")
        with open("Pele_ranking.txt", "w") as fout:
            fout.write("".join(output))
