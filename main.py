import os
import sys
sys.path.insert(0, "/sNow/easybuild/centos/7.4.1708/Skylake/software/schrodinger2017-4/internal/lib/python2.7/site-packages/")
import subprocess
import argparse
"""
NEW_LIB = '/sNow/easybuild/centos/7.4.1708/Skylake/software/schrodinger2017-4/internal/lib/ssl:/sNow/easybuild/centos/7.4.1708/Skylake/software/schrodinger2017-4/internal/lib:/sNow/easybuild/centos/7.4.1708/Skylake/software/schrodinger2017-4/mmshare-v4.0/lib/Linux-x86_64:/sNow/easybuild/centos/7.4.1708/Skylake/software/schrodinger2017-4/internal/lib/python2.7/site-packages'
try:
  OLD_LIB = os.environ["LD_LIBRARY_PATH"]
except KeyError:
  os.environ["LD_LIBRARY_PATH"] = NEW_LIB
else:
  os.environ["LD_LIBRARY_PATH"] = ':'.join([NEW_LIB, OLD_LIB])
"""
print(os.environ["LD_LIBRARY_PATH"])


import PlopRotTemp.main as plop
import Helpers.prepare_ligand as pl
import Helpers.check_env_var as env
import Helpers.pele_env as pele
import Adaptive.adaptive as ad
import Adaptive.clusterAdaptiveRun as cl
import Helpers.center_of_mass as cm 
import Helpers.constraints as ct
import SystemBuilder.system_prep as sp
import ppp.mutations_program as ppp

LIG_RES="LIG"
LIG_CHAIN="Z"
FORCEFIELD="OPLS2005"
PELE_CONFILE = "pele.conf"
CPUS = 3
RETURN = True
CLUSTERS = 40


ADAPTIVE_KEYWORDS = ["RESTART", "OUTPUT", "INPUT", "CPUS", "PELE_CFILE", "LIG_RES"]

PELE_KEYWORDS = ["NATIVE", "FORCEFIELD", "CHAIN", "CONSTRAINTS"]

NATIVE ='''
                        {{

                           "type": "rmsd",
                           
                           "Native": {{\n\
                            "path":\n\
                            "{}" }},\n\

                           "selection": {{ "chains": {{ "names": [ "$CHAIN" ] }} }},\n\

                           "includeHydrogens": false,\n\

                           "doSuperposition": false,\n\

                           "tag" : "ligandRMSD"\n\

                        }},\n\


'''

ADAPTIVE_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "Adaptive/clusterAdaptiveRun.py"))




def run(system, residue, chain, ligands, forcefield, confile, native, cpus, core, mtor, n, mae_charges, clean, only_plop):
    
    #Path Variables    
    pele_dir = os.path.abspath("{}_Pele".format(os.path.splitext(os.path.basename(system))[0]))
    adap_sh_input = os.path.join(pele_dir, "complex.pdb")
    adap_sh_output = os.path.join(pele_dir, "output_adaptive_short")
    pele_confile = os.path.join(pele_dir, PELE_CONFILE)
    cluster_output = os.path.join(pele_dir, "output_clustering")
    adap_l_input = "{}/initial_*"
    adap_l_output = os.path.join(pele_dir, "output_adaptive_long")

    #Template Variable
    native = NATIVE.format(os.path.abspath(native), chain) if native else native


    #Preparative for Pele
    receptor, ligand_mae, ligand_pdb = pl.prepare_ligand(system, residue, chain)
    cm_x, cm_y, cm_z = cm.center_of_mass(ligand_pdb)
    protein_constraints = ct.retrieve_constraints(system)
    builder = sp.SystemBuilder(ligands, receptor)
    structures = builder.structs()
    complexes = builder.systems()

    for structure, system in zip(structures, complexes):


	 system_fix, missing_residues = ppp.main(system)

 	 for residue, chain in missing_residues:
		
		 receptor, ligand_mae, ligand_pdb = pl.prepare_ligand(system, residue, chain)

	         template, rotamers_file = plop.main(ligand_mae, mtor, n, core, mae_charges, clean)
        
         ad_sh_temp, pele_temp, ad_l_temp = pele.set_pele_env(system, forcefield, template, rotamers_file, pele_dir)
 	
         ad.AdaptiveBuilder(pele_confile, PELE_KEYWORDS, native, forcefield, chain, "\n".join(protein_constraints))

       	 adaptive_short = ad.AdaptiveBuilder(ad_sh_temp, ADAPTIVE_KEYWORDS, RETURN, adap_sh_output, adap_sh_input, cpus, pele_confile, residue)

    	 adaptive_short.run()

         with pele.cd(adap_sh_output):

		cl.main(num_clusters=CLUSTERS, output_folder=cluster_output, ligand_resname=residue, atom_ids="")

       	 #adaptive_long = AdaptiveBuilder(ad_l_temp, ADAPTIVE_KEYWORDS, RETURN, adap_l_output, adap_l_input, cpus, pele_confile, residue)

      	 #adaptive_long.run()





if __name__ == "__main__":

    env.check_dependencies()
    
    parser = argparse.ArgumentParser(description='Run Adaptive Pele Platform')
    parser.add_argument('input', type=str, help='complex to run pele on')
    parser.add_argument('residue',  type=str, help='residue of the ligand to extract', default=LIG_RES)
    parser.add_argument('chain',  type=str, help='forcefield to use', default=LIG_CHAIN)
    parser.add_argument('ligands', type=str, help='ligands to run pele on')
    parser.add_argument('--forc',  type=str, help='chain of the ligand to extract', default=FORCEFIELD)
    parser.add_argument('--confile', type=str, help='your own pele configuration file', default=PELE_CONFILE)
    parser.add_argument('--native', type=str, help='native file to compare RMSD to', default="")
    parser.add_argument('--cpus', type=int, help='number of processors', default=CPUS)
    parser.add_argument("--core", type=int, help="Give one atom of the core section", default=-1)
    parser.add_argument("--mtor", type=int, help="Gives the maximum number of torsions allowed in each group.  Will freeze bonds to extend the core if necessary.", default=4)
    parser.add_argument("--n", type=int, help="Maximum Number of Entries in Rotamer File", default=1000)
    parser.add_argument("--mae_charges", help="Use charges in mae", action='store_true')
    parser.add_argument("--clean", help="Whether to clean up all the intermediate files", action='store_true')
    parser.add_argument("--only_plop", help="Whether to run PlopRotTemp or both", action='store_true')
    args = parser.parse_args()

    run(args.input, args.residue, args.chain, args.ligands, args.forc, args.confile, args.native, args.cpus, args.core, args.mtor, args.n, args.mae_charges, args.clean, args.only_plop)
