import os
import sys
import argparse
from PlopRotTemp.main import PlopRotTemp
from Helpers.prepare_ligand import prepare_ligand
from Helpers.check_env_var import check_dependencies
from Helpers.pele_env import set_pele_env
from Adaptive.adaptive import AdaptiveBuilder
#from Adaptive.clusterAdaptiveRun import cluster
# from Adaptive.long_adaptive import long_adaptive


LIG_RES="LIG"
LIG_CHAIN="Z"
FORCEFIELD="OPLS2005"
PELE_CONFILE = "pele.conf"
CPUS = 3
RETURN = True
CLUSTERS = 40


ADAPTIVE_KEYWORDS = ["RESTART", "OUTPUT", "INPUT", "CPUS", "PELE_CFILE", "LIG_RES"]

PELE_KEYWORDS = ["NATIVE", "FORCEFIELD", "CHAIN"]

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






def run(system, residue, chain, forcefield, confile, native, cpus, plop_opt, only_plop):

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


    #PlopRotTemp
    ligand_mae = prepare_ligand(system, residue, chain)
    template, rotamers_file = PlopRotTemp(ligand_mae, * plop_opt)
    
    #Pele
    if not only_plop:
        
        ad_sh_temp, pele_temp, ad_l_temp = set_pele_env(system, forcefield, template, rotamers_file, pele_dir)

        AdaptiveBuilder(pele_confile, PELE_KEYWORDS, native, forcefield, chain)

        adaptive_short = AdaptiveBuilder(ad_sh_temp, ADAPTIVE_KEYWORDS, RETURN, adap_sh_output, adap_sh_input, cpus, pele_confile, residue)

        adaptive_short.run()

        
        #cluster(num_clusters=CLUSTERS, outoutput_folder=cluster_output, ligand_resname=ligand_residue, inputFolder=adap_output)

        adaptive_long = AdaptiveBuilder(ad_l_temp, ADAPTIVE_KEYWORDS, RETURN, adap_l_output, adap_l_input, cpus, pele_confile, residue)

        adaptive_long.run()





if __name__ == "__main__":

    check_dependencies()
    
    parser = argparse.ArgumentParser(description='Run Adaptive Pele Platform')
    parser.add_argument('input', type=str, help='complex to run pele on')
    parser.add_argument('residue',  type=str, help='residue of the ligand to extract', default=LIG_RES)
    parser.add_argument('chain',  type=str, help='forcefield to use', default=LIG_CHAIN)
    parser.add_argument('--forc',  type=str, help='chain of the ligand to extract', default=FORCEFIELD)
    parser.add_argument('--confile', type=str, help='your own pele configuration file', default=PELE_CONFILE)
    parser.add_argument('--native', type=str, help='native file to compare RMSD to', default="")
    parser.add_argument('--cpus', type=int, help='number of processors', default=CPUS)
    parser.add_argument('--plop',  nargs='+', help='PlopRotTemp Options')
    parser.add_argument("--only_plop", help="Whether to run PlopRotTemp or both", action='store_true')
    args = parser.parse_args()

    run(args.input, args.residue, args.chain, args.forc, args.confile, args.native, args.cpus, args.plop, args.only_plop)