import matplotlib
matplotlib.use('Agg')
from pyemma import config
config.show_progress_bars = False
#config.mute = True
import os
import msm_pele.Helpers.check_env_var as env
env.check_dependencies()
import shutil
import random
import glob
import argparse
import msm_pele.constants as cs
import msm_pele.PlopRotTemp.launcher as plop
import msm_pele.Helpers.helpers as hp
import msm_pele.Helpers.pele_env as pele
import msm_pele.Helpers.center_of_mass as cm
import msm_pele.Helpers.simulation as ad
import msm_pele.Helpers.clusterAdaptiveRun as cl
import msm_pele.Helpers.system_prep as sp
import msm_pele.Box.box as bx
import msm_pele.PPP.mut_prep4pele as ppp
import msm_pele.Analysis.msm_analysis as msm
import msm_pele.Helpers.missing_residues as mr
import msm_pele.Helpers.solventOBCParamsGenerator as obc

__version__ = "2.0.0"

def run(args):
    # Build folders and logging
    env = pele.EnviroBuilder.build_env(args)

    if args.restart == "all":

        # Build System
        env.logger.info("Checking {} system for Pele".format(args.residue))
        syst = sp.SystemBuilder.build_system(args.system, args.mae_lig, args.residue, env.pele_dir)

        # Prepare System
        system_fix, missing_residues, gaps, metals, protein_constraints = ppp.main(syst.system, env.pele_dir, charge_terminals=args.charge_ter,
                no_gaps_ter=args.gaps_ter, mid_chain_nonstd_residue=env.nonstandard, renumber=env.renumber, dynamic_waters=env.water)
        env.logger.info(cs.SYSTEM.format(system_fix, missing_residues, gaps, metals))

        # Parametrize Ligand
        env.logger.info("Creating template for residue {}".format(args.residue))
        with hp.cd(env.pele_dir):
        	plop.parametrize_miss_residues(args, env, syst)
        env.logger.info("Template {}z created".format(args.residue.lower()))

        # Parametrize missing residues
        for res, __, _ in missing_residues:
            if res != args.residue:
                env.logger.info("Creating template for residue {}".format(res))
                with hp.cd(env.pele_dir):
                    mr.create_template(args, env)
                    env.logger.info("Template {}z created".format(res))

        # Parametrize solvent parameters if need it
        if env.solvent == "OBC":
            shutil.copy(env.obc_tmp, env.obc_file)
            for template in glob.glob(os.path.join(env.template_folder, "*")):
                obc.main(template, env.obc_file)

        # Fill in Simulation Templates
        ad.SimulationBuilder(env.pele_exit_temp,  env.topology, cs.EX_PELE_KEYWORDS, env.native, args.forcefield, args.chain, "\n".join(protein_constraints), env.license, env.log, env.solvent, env.dynamic_water, env.temp)
        ad.SimulationBuilder(env.pele_temp,  env.topology, cs.EX_PELE_KEYWORDS, env.native, args.forcefield, args.chain, "\n".join(protein_constraints),  env.license, env.log, env.solvent, env.dynamic_water, env.temp)



    if args.restart in ["all", "pele"]:

        #Check restart variables
        try:
            initial_iteration = max([ int(os.path.basename(os.path.normpath(folder))) for folder in glob.glob(os.path.join(env.adap_l_output, "*/")) if os.path.basename(os.path.normpath(folder)).isdigit() ])
        except ValueError: 
            initial_iteration = 0

        # For each iteration
        for i in range(initial_iteration, env.iterations):

            if not env.one_exit or i == 0:
                # Run Multiple Adaptive Exit if no specified the contrary
                env.logger.info("Running ExitPath Adaptive {}".format(i))
                env.update_variable_for_iteration(i)
                shutil.copy(env.adap_exit_template, env.ad_ex_temp)
                simulation = ad.SimulationBuilder(env.ad_ex_temp, env.topology, cs.EX_ADAPTIVE_KEYWORDS, cs.RESTART, env.adap_ex_output, env.adap_ex_input, env.cpus, env.pele_exit_temp, env.residue, env.equil_steps, env.random_num, env.exit_iters, env.eq_struct, env.usesrun)
                simulation.run_adaptive(env, hook=True)
                env.logger.info("ExitPath Adaptive run successfully")
                

                # KMeans Clustering with different seed
                env.logger.info("Running Exit Path Clustering")
                with hp.cd(env.adap_ex_output):
                    seed = 742304 if env.test else random.randint(1,1000000)
                    env.logger.info(seed)
                    cluster_centers = cl.main(env.clusters, env.cluster_output, args.residue, "", env.cpus, env.adap_ex_output,
                       env.topology, env.sasamin, env.sasamax, env.sasa, env.perc_sasa_min, env.perc_sasa_int, seed=seed, iteration=i)
                env.logger.info("Exit Path Clustering run successfully")

            # Create Exploration Box
            env.logger.info("Creating box")
            box, BS_sasa_min, BS_sasa_max = bx.create_box(cluster_centers, env, i)
            env.logger.info("Box created successfully ")

            # Explore
            env.logger.info("Running standard Pele")
            inputs = [ cs.INPUT_PELE.format(f) for f in glob.glob(env.adap_l_input) ]
            output = os.path.join(env.adap_l_output, str(i))
            if not os.path.isdir(output):
                os.mkdir(output)
            hp.change_output(env.pele_temp, output)
            simulation = ad.SimulationBuilder(env.pele_temp,  env.topology, cs.PELE_KEYWORDS, cs.RESTART, os.path.join(env.adap_l_output, output),
                ",\n".join(inputs), env.random_num, env.steps, box, env.box_metric, BS_sasa_min, BS_sasa_max, env.temp)
            time_sim = simulation.run_pele(env, limitTime=env.time)
            env.logger.info("Pele run successfully in {}".format(time_sim))

            # Analyse simulation
            env.logger.info("Running MSM analysis")
            msm.analyse_results(env, runTica=False, last=i == (env.iterations-1))
            env.logger.info("MSM analysis run successfully")

    if args.restart in ["msm", "analyse"]:

        # MSM Final Analysis
        env.logger.info("Running MSM analysis")
        env.update_variable_for_iteration(env.iterations-1)
        msm.analyse_results(env, runTica=False, last=True)
        env.logger.info("MSM analysis run successfully")

    env.logger.info("{} System run successfully".format(args.residue))


def parse_args(args=[]):

    parser = argparse.ArgumentParser(description='Run Adaptive Pele Platform')
    parser.add_argument('system', type=str, help='complex to run pele on')
    parser.add_argument('residue', type=str, help='residue of the ligand to extract', default=cs.LIG_RES)
    parser.add_argument('chain', type=str, help='chain of the ligand to extract', default=cs.LIG_CHAIN)
    parser.add_argument("--mae_lig", type=str, help="ligand .mae file to include QM charges coming from jaguar")
    parser.add_argument("--box", type=str, help="Exploration box for Pele")
    parser.add_argument("--charge_ter", help="Charge protein terminals", action='store_true')
    parser.add_argument("--gaps_ter", help="Include TER when a possible gap is found", action='store_true')
    parser.add_argument("--clust", type=int, help="Numbers of clusters to start PELE's exploration with", default=cs.CLUSTERS)
    parser.add_argument('--forcefield', '-f', type=str, help='chain of the ligand to extract', default=cs.FORCEFIELD)
    parser.add_argument('--confile', type=str, help='your own pele configuration file', default=cs.PELE_CONFILE)
    parser.add_argument('--native', type=str, help='native file to compare RMSD to', default="")
    parser.add_argument('--cpus', type=int, help='number of processors', default=cs.CPUS)
    parser.add_argument("--core", type=int, help="Give one atom of the core section", default=-1)
    parser.add_argument("--mtor", type=int, help="Gives the maximum number of torsions allowed in each group.  Will freeze bonds to extend the core if necessary.", default=4)
    parser.add_argument("--n", type=int, help="Maximum Number of Entries in Rotamer File", default=1000)
    parser.add_argument("--clean", help="Whether to clean up all the intermediate files", action='store_true')
    parser.add_argument("--restart", type=str, help="Restart the platform from [all, adaptive, pele, msm, analyse] with these keywords", default=cs.PLATFORM_RESTART)
    parser.add_argument("--gridres", type=str, help="Rotamers angle resolution", default=cs.GRIDRES)
    parser.add_argument("--precision", action='store_true', help="Use a more agressive control file to achieve better convergence")
    parser.add_argument("--precision2", action='store_true', help="Use an intermediate control file to achieve better convergence")
    parser.add_argument("--test", action='store_true', help="Run a fast msm_pele test")
    parser.add_argument("--user_center", "-c", nargs='+', type=float, help='center of the box', default=None)
    parser.add_argument("--user_radius", "-r", nargs='+', type=float,  help="Radius of the box", default=None)
    parser.add_argument("--folder", "-wf", type=str,  help="Folder to apply the restart to", default=None)
    parser.add_argument("--pdb", action='store_true',  help="Use pdb files as output")
    parser.add_argument("--nonstandard", nargs="+",  help="Mid Chain non standard residues to be treated as ATOM not HETATOM", default = [])
    parser.add_argument("--lagtime", type=int,  help="MSM Lagtime to use", default=100)
    parser.add_argument("--steps", type=int,  help="MSM Steps to use", default=1000)
    parser.add_argument("--temp", type=int,  help="Temperature to use", default=1000)
    parser.add_argument("--msm_clust", type=int,  help="Number of clusters created to converge MSM", default=200)
    parser.add_argument("--time", type=int,  help="Limit of time to run pele exploration", default=None)
    parser.add_argument("--log", action='store_true',  help="Print LogFiles when running PELE")
    parser.add_argument("--nonrenum", action='store_false',  help="Don't renumber structure")
    parser.add_argument("--sasa", nargs="+",  type=float, help="Interval of sasa used for adaptive exit clusterization", default=[])
    parser.add_argument("--box_type",  type=str, help="Type of box to use. [multiple (default), fixed]", default="multiple")
    parser.add_argument("--ext_temp",  nargs="+", type=str, help="Use external template to parametrize the ligand i.e. /path/mgz", default=[])
    parser.add_argument("--nosasa",  action='store_true', help="Do not filter clusters by sasa i.e. --nosasa")
    parser.add_argument("--perc_sasa",  nargs="+", type=float, help="Distribution of clusters at the adaptive exit. default [0.25, 0.5, 0.25] i.e. 0.1 0.7 0.1", default = [0.25, 0.5, 0.25])
    parser.add_argument("--box_metric",  action='store_true', help="Do not filter clusters by sasa i.e. --nosasa")
    parser.add_argument('--iterations', type=int, help='Number of MSM iterations. default [1] i.e. --iterations 3', default=1)
    parser.add_argument('--exit_iters', type=int, help='Number of exit epochs on simulation. default [100] i.e. --exit_iters 1', default=100)
    parser.add_argument('--eq_struct', type=int, help='Number of structures extracted from inital equilibration to launch exit simulation. default [10] i.e. --eq_struct 1', default=10)
    parser.add_argument('--solvent', type=str, help='Type of implicit solvent (OBC/VDGBNP). default [OBC]. i.e. --solvent VDGBNP', default="OBC")
    parser.add_argument("--one_exit",  action='store_true', help="Perform only one adaptive exit simulation")
    parser.add_argument("--noRMSD",  action='store_true', help="De not use keep track RMSD to the initial position. NOT TESTED.")
    parser.add_argument("--water",  nargs="+", help="Waters to perturb. i.e B:1 B:2 B:3")
    parser.add_argument("--water_radius",  type=float, help="Radius of the box to define the exploration area of MC waters. i.e --water_radius 7 default=10", default=10)
    parser.add_argument("--water_temp",  type=int, help="Temperature of water MC. i.e --water_temp 1000 default=500", default=500)
    parser.add_argument("--water_constr",  type=float, help="Constraint on the waters MC. i.e ----water_const 0.5 default=0.2", default=0.2)
    parser.add_argument("--water_trials",  type=int, help="Steric trials on the waters MC. i.e --water_trials 2000 default=1000", default=1000)
    
    args = parser.parse_args(args) if args else parser.parse_args()
    return args

if __name__ == "__main__":
    args = parse_args()
    if(args.clust > args.cpus and args.restart != "msm" and not args.test ):
        raise ValueError(cs.CLUSTER_ERROR.format(args.cpus, args.clust))
    else:
        run(args)
