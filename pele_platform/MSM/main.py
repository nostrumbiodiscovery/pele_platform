import os
import pele_platform.Utilities.Helpers.check_env_var as env
env.check_dependencies()
import argparse
import pele_platform.constants as cs
import pele_platform.Utilities.PlopRotTemp.launcher as plop
import pele_platform.Utilities.Helpers.helpers as hp
import pele_platform.Utilities.Parameters.pele_env as pele
import pele_platform.Utilities.Helpers.simulation as ad
import pele_platform.MSM.clusterAdaptiveRun as cl
import pele_platform.Utilities.Helpers.system_prep as sp
import pele_platform.MSM.box as bx
import pele_platform.Utilities.PPP.main as ppp
import pele_platform.MSM.msm_analysis as msm
import pele_platform.Utilities.Helpers.missing_residues as mr

__version__ = "1.0.3"

def run(args):
    # Build folders and logging
    env = pele.EnviroBuilder.build_env(args)

    if args.restart in cs.FIRST_RESTART:

        # Build System
        env.logger.info("Checking {} system for Pele".format(args.residue))
        syst = sp.SystemBuilder.build_system(args.system, args.mae_lig, args.residue, env.pele_dir)

        # Prepare System
        system_fix, missing_residues, gaps, metals, env.constraints = ppp.main(syst.system, env.pele_dir, charge_terminals=args.charge_ter, no_gaps_ter=args.gaps_ter)
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

    if args.restart in cs.SECOND_RESTART:
        # Run Adaptive Exit
        env.logger.info("Running ExitPath Adaptive")
        adaptive_exit = ad.SimulationBuilder(env.ad_ex_temp, env.pele_exit_temp, env)
        adaptive_exit.run(hook=True)
        env.logger.info("ExitPath Adaptive run successfully")


    if args.restart in cs.THIRD_RESTART:

        #KMeans Clustering
        if not os.path.isfile(env.clusters_output):
            env.logger.info("Running MSM Clustering")
            with hp.cd(env.adap_ex_output):
                cl.main(env.clusters, env.cluster_output, args.residue, "", env.cpus, env.topology)
            env.logger.info("MSM Clustering run successfully")
        else:
            pass

        # Create Box
        env.logger.info("Creating box")
        env.box_center, env.box_radius, env.sasa_min, env.sasa_max = bx.create_box(args, env)
        env.logger.info("Box with center {} radius {} was created".format(env.box_center, env.box_radius))

        # Pele Exploration
        env.logger.info("Running standard Pele")
        adaptive_long = ad.SimulationBuilder(env.ad_l_temp, env.pele_temp, env)
        adaptive_long.run()
        env.logger.info("Pele run successfully")

    if args.restart in cs.FOURTH_RESTART:

        # MSM Analysis
        env.logger.info("Running MSM analysis")
        msm.analyse_results(env, args, runTica=False)
        env.logger.info("MSM analysis run successfully")

        env.logger.info("{} System run successfully".format(args.residue))

    return env

if __name__ == "__main__":

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
    parser.add_argument("--restart", type=str, help="Restart the platform from [all, pele, msm] with these keywords", default=cs.PLATFORM_RESTART)
    parser.add_argument("--gridres", type=str, help="Rotamers angle resolution", default=cs.GRIDRES)
    parser.add_argument("--precision", action='store_true', help="Use a more agressive control file to achieve better convergence")
    parser.add_argument("--test", action='store_true', help="Run a fast pele_platform test")
    parser.add_argument("--user_center", "-c", nargs='+', type=float, help='center of the box', default=None)
    parser.add_argument("--user_radius", "-r", type=float,  help="Radius of the box", default=None)
    parser.add_argument("--folder", "-wf", type=str,  help="Folder to apply the restart to", default=None)
    parser.add_argument("--pdb", action='store_true',  help="Use pdb files as output")
    
    args = parser.parse_args()
    if(args.clust > args.cpus and args.restart != "msm" and not args.test ):
        raise ValueError(cs.CLUSTER_ERROR.format(args.cpus, args.clust))
    else:
        run(args)
