import os
import MSM_PELE.Helpers.check_env_var as env
env.check_dependencies()
import logging
import argparse
import random
import MSM_PELE.constants as cs
import MSM_PELE.PlopRotTemp.main as plop
import MSM_PELE.Helpers.helpers as hp
import MSM_PELE.Helpers.pele_env as pele
import MSM_PELE.Helpers.simulation as ad
import MSM_PELE.Helpers.clusterAdaptiveRun as cl
import MSM_PELE.Helpers.center_of_mass as cm
import MSM_PELE.Helpers.constraints as ct
import MSM_PELE.Helpers.system_prep as sp
import MSM_PELE.Helpers.box as bx
import MSM_PELE.PPP.mut_prep4pele as ppp
import MSM_PELE.Helpers.msm_analysis as msm

# DEFAULT VALUES
COMPLEX = "complex.pdb"
RESULTS = "results"
LIG_RES = "LIG"
LIG_CHAIN = "Z"
FORCEFIELD = "OPLS2005"
PELE_CONFILE = "pele.conf"
CPUS = 140
RESTART = True
CLUSTERS = 40
PLATFORM_RESTART = "all"
EQ_STEPS = 750

# KEYWORDS
ADAPTIVE_KEYWORDS = ["RESTART", "OUTPUT", "INPUT", "CPUS", "PELE_CFILE", "LIG_RES", "SEED"]
EX_ADAPTIVE_KEYWORDS = ["RESTART", "OUTPUT", "INPUT", "CPUS", "PELE_CFILE", "LIG_RES", "EQ_STEPS", "SEED"]
EX_PELE_KEYWORDS = ["NATIVE", "FORCEFIELD", "CHAIN", "CONSTRAINTS", "CPUS", "LICENSES"]
PELE_KEYWORDS = ["BOX_CENTER", "BOX_RADIUS"]
NATIVE = '''
                        {{

                           "type": "rmsd",

                           "Native": {{\n\
                            "path":\n\
                            "{}" }},\n\

                           "selection": {{ "chains": {{ "names": [ "{}" ] }} }},\n\

                           "includeHydrogens": false,\n\

                           "doSuperposition": false,\n\

                           "tag" : "ligandRMSD"\n\

                        }},\n\


'''

SYSTEM = "System {} checked successfully\n\t**Missing residues found {}\n\t**Gaps found {}\n\t**Metals found {}"



# FOLDERS&PATH
DIR = os.path.dirname(__file__)
ADAPTIVE_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "Adaptive/clusterAdaptiveRun.py"))
FOLDERS = ["",
           "DataLocal/Templates/OPLS2005/HeteroAtoms/",
           "DataLocal/Templates/AMBER99sb/HeteroAtoms/",
           "DataLocal/Templates/AMBER99sbBSC0/HeteroAtoms/",
           "DataLocal/LigandRotamerLibs"]

# ERRORS
CLUSTER_ERROR = "Number of cpus ({}) must be bigger than clusters ({})"





def run(system, residue, chain, mae_lig, charge_ter, gaps_ter, clusters, forcefield, confile, native, cpus, core, mtor, n, clean, restart):

    # Preparative for Pele
    template = None
    rotamers_file = None
    random_num = random.randrange(1,70000) 
    license = '''"{}"'''.format(cs.LICENSE)
    equil_steps = int(EQ_STEPS/cpus)
    pele_dir = os.path.abspath("{}_Pele".format(residue))

    if clusters > cpus:
        raise ValueError(CLUSTER_ERROR.format(cpus, clusters))

    if restart == "all":
        pele_dir = pele.is_repited(pele_dir)
    else:
        pele_dir = pele.is_last(pele_dir)
    logger, log_name = hp.set_logger(pele_dir, residue)
    native = NATIVE.format(os.path.abspath(native), chain) if native else native
    system_fix = "{}_processed.pdb".format(os.path.abspath(os.path.splitext(system)[0]))
    adap_ex_input = os.path.join(pele_dir, os.path.basename(system_fix))
    adap_ex_output = os.path.join(pele_dir, "output_adaptive_exit")
    cluster_output = os.path.join(pele_dir, "output_clustering")
    adap_l_input = "{}/initial_*"
    adap_l_output = os.path.join(pele_dir, "output_pele")
    ad_ex_temp = os.path.join(pele_dir, "adaptive_exit.conf")
    ad_l_temp = os.path.join(pele_dir, "adaptive_long.conf")
    pele_exit_temp = os.path.join(pele_dir, "pele_exit.conf")
    pele_temp = os.path.join(pele_dir, "pele.conf")
    box_temp = os.path.join(pele_dir, "box.pdb")
    clusters_output = os.path.join(cluster_output, "clusters_40_KMeans_allSnapshots.pdb")
    lig_ref = os.path.join(pele_dir, "ligand.pdb")

    if restart == "all":

        # Building system and ligand
        logger.info("Checking {} system for Pele".format(residue))
        if mae_lig:
            receptor = system
            lig = mae_lig
            lig_ref = sp.convert_pdb(mae_lig)
            system = sp.build_complex(receptor, lig_ref)
        else:
            receptor, lig_ref = sp.retrieve_receptor(system, residue)
            lig, residue = sp.convert_mae(lig_ref)
        system_fix, missing_residues, gaps, metals = ppp.main(system, charge_terminals=charge_ter, no_gaps_ter=gaps_ter)
        protein_constraints = ct.retrieve_constraints(system_fix, gaps, metals)
        logger.info(SYSTEM.format(system_fix,missing_residues, gaps, metals))

        for res, _, _ in missing_residues:
            logger.info("Creating template for residue {}".format(res))
            if mae_lig:
                mae_charges = True
                template, rotamers_file = plop.main(mae_lig, mtor, n, core, mae_charges, clean)
                hp.silentremove([system])
            else:
                mae_charges = False
                template, rotamers_file = plop.main(lig, mtor, n, core, mae_charges, clean)
                hp.silentremove([lig])
            logger.info("Template {} created".format(template))

        files_to_copy = [os.path.join(DIR, "Templates/box.pdb"), os.path.join(DIR, "Templates/pele.conf"),
                 os.path.join(DIR, "Templates/adaptive_exit.conf"), os.path.join(DIR, "Templates/adaptive_long.conf"),
                 os.path.join(DIR, "Templates/pele_exit.conf")]
        files_to_move = [system_fix, log_name, lig_ref]
        directories = FOLDERS
        directories.extend(["output_pele", "output_adaptive_exit", "output_clustering"]) 
        pele.set_pele_env(directories, files_to_copy, forcefield, template, rotamers_file, files_to_move, pele_dir)
        ad.SimulationBuilder(pele_exit_temp, EX_PELE_KEYWORDS, native, forcefield, chain, "\n".join(protein_constraints), cpus, license)
        ad.SimulationBuilder(pele_temp, EX_PELE_KEYWORDS, native, forcefield, chain, "\n".join(protein_constraints), cpus, license)

        logger.info("Running ExitPath Adaptive")
        adaptive_exit = ad.SimulationBuilder(ad_ex_temp, EX_ADAPTIVE_KEYWORDS, RESTART, adap_ex_output, adap_ex_input, cpus, pele_exit_temp, residue, equil_steps, random_num)
        adaptive_exit.run()
        logger.info("ExitPath Adaptive run successfully")


    if restart in ["all", "pele"]:

        if not clusters_output:
            logger.info("Running MSM Clustering")
            with hp.cd(adap_ex_output):
                cl.main(num_clusters=clusters, output_folder=cluster_output, ligand_resname=residue, atom_ids="")
            logger.info("MSM Clustering run successfully")
        else:
            pass
 
        logger.info("Creating box")
        center_mass = cm.center_of_mass(lig_ref)
        center, radius = bx.main(adap_ex_output, clusters_output, center_mass)
        bx.build_box(center, radius, box_temp)
        logger.info("Box Created")

        logger.info("Running standard Pele")
        ad.SimulationBuilder(pele_temp, PELE_KEYWORDS, center, radius)
        adaptive_long = ad.SimulationBuilder(ad_l_temp, ADAPTIVE_KEYWORDS, RESTART, adap_l_output, adap_l_input, cpus, pele_temp, residue, random_num)
        adaptive_long.run()
        logger.info("Pele run successfully")

    if restart in ["all", "pele", "msm"]:

        logger.info("Running MSM analysis")
        msm.analyse_results(adap_l_output, residue, cpus, pele_dir)
        logger.info("MSM analysis run successfully")        

        logger.info("{} System run successfully".format(residue))


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Run Adaptive Pele Platform')
    parser.add_argument('input', type=str, help='complex to run pele on')
    parser.add_argument('residue', type=str, help='residue of the ligand to extract', default=LIG_RES)
    parser.add_argument('chain', type=str, help='chain of the ligand to extract', default=LIG_CHAIN)
    parser.add_argument("--mae_lig", type=str, help="ligand .mae file to include QM charges coming from jaguar")
    parser.add_argument("--charge_ter", help="Charge protein terminals", action='store_true')
    parser.add_argument("--gaps_ter", help="Include TER when a possible gap is found", action='store_true')
    parser.add_argument("--clust", type=int, help="Numbers of clusters to start PELE's exploration with", default=CLUSTERS)
    parser.add_argument('--forc', type=str, help='chain of the ligand to extract', default=FORCEFIELD)
    parser.add_argument('--confile', type=str, help='your own pele configuration file', default=PELE_CONFILE)
    parser.add_argument('--native', type=str, help='native file to compare RMSD to', default="")
    parser.add_argument('--cpus', type=int, help='number of processors', default=CPUS)
    parser.add_argument("--core", type=int, help="Give one atom of the core section", default=-1)
    parser.add_argument("--mtor", type=int, help="Gives the maximum number of torsions allowed in each group.  Will freeze bonds to extend the core if necessary.", default=4)
    parser.add_argument("--n", type=int, help="Maximum Number of Entries in Rotamer File", default=1000)
    parser.add_argument("--clean", help="Whether to clean up all the intermediate files", action='store_true')
    parser.add_argument("--restart", type=str, help="Restart the platform from [all, pele, msm] with these keywords", default=PLATFORM_RESTART)
    args = parser.parse_args()

    if args.clust > args.cpus:
        raise ValueError(CLUSTER_ERROR.format(args.cpus, args.clust))
    else:
        run(args.input, args.residue, args.chain, args.mae_lig, args.charge_ter, args.gaps_ter, args.clust, args.forc, args.confile, args.native, args.cpus, args.core, args.mtor, args.n, args.clean, args.restart)
