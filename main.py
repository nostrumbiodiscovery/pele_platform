import sys
import os
import logging
import argparse
import PlopRotTemp.main as plop
import Helpers.helpers as hp
import Helpers.prepare_ligand as pl
import Helpers.check_env_var as env
import Helpers.pele_env as pele
import Adaptive.adaptive as ad
import Adaptive.clusterAdaptiveRun as cl
import Helpers.center_of_mass as cm
import Helpers.constraints as ct
import SystemBuilder.system_prep as sp
import SystemBuilder.box as bx
import ppp.mutations_program as ppp
import msm.analysis as msm

sys.path.insert(0, "/sNow/easybuild/centos/7.4.1708/Skylake/software/schrodinger2017-4/internal/lib/python2.7/site-packages/")


COMPLEX = "complex.pdb"
RESULTS = "results"
LIG_RES = "LIG"
LIG_CHAIN = "Z"
FORCEFIELD = "OPLS2005"
PELE_CONFILE = "pele.conf"
CPUS = 3
RESTART = True
CLUSTERS = 40

ADAPTIVE_KEYWORDS = ["RESTART", "OUTPUT", "INPUT", "CPUS", "PELE_CFILE", "LIG_RES"]

EX_PELE_KEYWORDS = ["NATIVE", "FORCEFIELD", "CHAIN", "CONSTRAINTS"]

PELE_KEYWORDS = ["NATIVE", "FORCEFIELD", "CHAIN", "CONSTRAINTS", "BOX_CENTER", "BOX_RADIUS"]

NATIVE = '''
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

# Output Constants
RANKING_FILE = "Pele_ranking.txt"
OUTPUT_HEADER = "#Residue Epoch DG StdDG Db StdDb Conv\n#==================================\n"

# Folders and files
FOLDERS = ["",
           "DataLocal/Templates/OPLS2005/HeteroAtoms/",
           "DataLocal/Templates/AMBER99sb/HeteroAtoms/",
           "DataLocal/Templates/AMBER99sbBSC0/HeteroAtoms/",
           "DataLocal/LigandRotamerLibs",
           ]


# Log Constants
LOG_FILENAME = "output.out"
LOG_FORMAT = "%(asctime)s:%(levelname)s:%(message)s"

# directory
DIR = os.path.dirname(__file__)

# Logging definition block
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
formatter = logging.Formatter(LOG_FORMAT)
file_handler = logging.FileHandler(LOG_FILENAME, mode='w')
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)


def run(system, residue, chain, ligands, forcefield, confile, native, cpus, core, mtor, n, mae_charges, clean, only_plop):

    # Template Variable
    native = NATIVE.format(os.path.abspath(native), chain) if native else native
    protein_constraints = ct.retrieve_constraints(system)
    r_cm = cm.center_of_mass(ligand_pdb)

    # Exit pathway files
    logger.info("Preparing ExitPath Adaptive Env")
    pele_dir = os.path.abspath("{}_Pele".format(os.path.splitext(system)[0]))
    pele_confile = os.path.join(pele_dir, PELE_CONFILE)
    adap_ex_input = os.path.join(pele_dir, os.path.basename(system))
    adap_ex_output = os.path.join(pele_dir, RESULTS)
    ad_ex_temp = os.path.join(pele_dir, "adaptive_exit.conf")
    # Running ppp
    system_fix, missing_residues = ppp.main(system)
    # Files to create
    ex_files = [system_fix, os.path.join(DIR, "PeleTemplates/pele.conf"), os.path.join(DIR, "PeleTemplates/adaptive_exit.conf")]
    ex_folders = FOLDERS
    ex_folders.append(RESULTS)
    # Look for exit pathway with Adaptive PELE
    receptor, ligand_mae, ligand_pdb = pl.prepare_ligand(system_fix, residue, chain)
    template, rotamers_file = plop.main(ligand_mae, mtor, n, core, mae_charges, clean)
    logger.info("Running exit Adaptive")
    pele.set_pele_env(system_fix, ex_folders, ex_files, forcefield, template, rotamers_file, pele_dir)
    ad.AdaptiveBuilder(pele_confile, EX_PELE_KEYWORDS, native, forcefield, chain, "\n".join(protein_constraints))
    adaptive_exit = ad.AdaptiveBuilder(ad_ex_temp, ADAPTIVE_KEYWORDS, RESTART, adap_ex_output, adap_ex_input, cpus, pele_confile, residue)
    adaptive_exit.run()

    # Preparative for Pele
    logger.info("Retrieving Ligands & Complexes")
    receptor = sp.retrieve_receptor(system)
    ligands, complexes, residues = sp.build_complexes(ligands, receptor)
    pele_dirs = [os.path.abspath("{}_Pele".format(residue)) for residue in residues]

    # Run PlopRotTemp + Pele
    for residue, lig, sys, pele_dir in zip(residues, ligands, complexes, pele_dirs):

        # Path Variables
        logger.info("Creating Pele env")
        adap_sh_input = os.path.join(pele_dir, COMPLEX)
        adap_sh_output = os.path.join(pele_dir, "output_adaptive_short")
        cluster_output = os.path.join(pele_dir, "output_clustering")
        adap_l_input = "{}/initial_*"
        adap_l_output = os.path.join(pele_dir, "output_adaptive_long")
        ad_sh_temp = os.path.join(pele_dir, "adaptive_short.conf")
        ad_l_temp = os.path.join(pele_dir, "adaptive_long.conf")
        pele_temp = os.path.join(pele_dir, "pele.conf")
        box_temp = os.path.join(pele_dir, "box.pdb")

        files = [COMPLEX, os.path.join(DIR, "PeleTemplates/box.pdb"), os.path.join(DIR, "PeleTemplates / pele.conf"),
                 os.path.join(DIR, "PeleTemplates/adaptive_short.conf"), os.path.join(DIR, "PeleTemplates/adaptive_long.conf")]
        directories = FOLDERS
        directories.extend(["output_adaptive_long", "output_adaptive_short", "output_clustering"])

        # Fix protein
        system_fix, missing_residues = ppp.main(sys)
        # Produce Templates of all missing residues
        for res, chain in missing_residues:
            logger.info("Creating template for residue {}".format(res))
            template, rotamers_file = plop.main(lig, mtor, n, core, mae_charges, clean)
            hp.silentremove(lig)

        # Run Pele
        logger.info("Setting Pele environment folders")
        pele.set_pele_env(system_fix, directories, files, forcefield, template, rotamers_file, pele_dir)
        logger.info("Create box")
        # HOWWW???
        box = bx.build_box(CENTER, RADIUS, box_temp)
        logger.info("Creating pele control file")
        ad.AdaptiveBuilder(pele_temp, PELE_KEYWORDS, native, forcefield, chain, "\n".join(protein_constraints), box.center, box.r)
        logger.info("Creating adaptive control file")
        adaptive_short = ad.AdaptiveBuilder(ad_sh_temp, ADAPTIVE_KEYWORDS, RESTART, adap_sh_output, adap_sh_input, cpus, pele_confile, residue)
        logger.info("Running Adaptive")
        adaptive_short.run()
        logger.info("MSM Clustering")
        with pele.cd(adap_sh_output):
            cl.main(num_clusters=CLUSTERS, output_folder=cluster_output, ligand_resname=residue, atom_ids="")
        logger.info("Running standard Pele")
        adaptive_long = ad.AdaptiveBuilder(ad_l_temp, ADAPTIVE_KEYWORDS, RESTART, adap_l_output, adap_l_input, cpus, pele_confile, residue)
        adaptive_long.run()
    msm.analyse_results(adap_l_output, residue)

    # Analyze results
    output = msm.summerize(pele_dirs, residues)
    output.insert(0, OUTPUT_HEADER)
    with open(RANKING_FILE, "w") as fout:
        fout.write("".join(output))


if __name__ == "__main__":

    env.check_dependencies()

    parser = argparse.ArgumentParser(description='Run Adaptive Pele Platform')
    parser.add_argument('input', type=str, help='complex to run pele on')
    parser.add_argument('residue', type=str, help='residue of the ligand to extract', default=LIG_RES)
    parser.add_argument('chain', type=str, help='forcefield to use', default=LIG_CHAIN)
    parser.add_argument('ligands', type=str, help='ligands to run pele on')
    parser.add_argument('--forc', type=str, help='chain of the ligand to extract', default=FORCEFIELD)
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
