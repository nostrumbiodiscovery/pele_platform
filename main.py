import os
import MSM_PELE.Helpers.check_env_var as env
env.check_dependencies()
import shutil
import argparse
import MSM_PELE.constants as cs
import MSM_PELE.PlopRotTemp.main as plop
import MSM_PELE.Helpers.helpers as hp
import MSM_PELE.Helpers.pele_env as pele
import MSM_PELE.Helpers.simulation as ad
import MSM_PELE.Helpers.clusterAdaptiveRun as cl
import MSM_PELE.Helpers.center_of_mass as cm
import MSM_PELE.Helpers.system_prep as sp
import MSM_PELE.Helpers.box as bx
import MSM_PELE.PPP.mut_prep4pele as ppp
import MSM_PELE.Helpers.msm_analysis as msm
import MSM_PELE.Helpers.missing_residues as mr

__version__ = "1.0.3"

def run(args):
    # Build folders and logging
    env = pele.EnviroBuilder.build_env(args)

    if args.restart == "all":

        # Build System
        env.logger.info("Checking {} system for Pele".format(args.residue))
        syst = sp.SystemBuilder.build_system(args.system, args.mae_lig, args.residue, env.pele_dir)

        # Prepare System
        system_fix, missing_residues, gaps, metals, protein_constraints = ppp.main(syst.system, env.pele_dir, charge_terminals=args.charge_ter, no_gaps_ter=args.gaps_ter)
        env.logger.info(cs.SYSTEM.format(system_fix, missing_residues, gaps, metals))

        # Parametrize Ligand
        env.logger.info("Creating template for residue {}".format(args.residue))
        plop.parametrize_miss_residues(args, env, syst)
        env.logger.info("Template {}z created".format(args.residue.lower()))

        # Parametrize missing residues
        for res, __, _ in missing_residues:
            if res != args.residue:
                env.logger.info("Creating template for residue {}".format(res))
                mr.create_template(system_fix, res, env.pele_dir, args.forcefield)
                env.logger.info("Template {}z created".format(res))

        # Fill in Simulation Templates
        ad.SimulationBuilder(env.pele_exit_temp,  env.topology, cs.EX_PELE_KEYWORDS, env.native, args.forcefield, args.chain, "\n".join(protein_constraints), args.cpus, env.license)
        ad.SimulationBuilder(env.pele_temp,  env.topology, cs.EX_PELE_KEYWORDS, env.native, args.forcefield, args.chain, "\n".join(protein_constraints), args.cpus, env.license)

    if args.restart in ["all", "adaptive"]:
        # Run Adaptive Exit
        env.logger.info("Running ExitPath Adaptive")
        adaptive_exit = ad.SimulationBuilder(env.ad_ex_temp, env.topology, cs.EX_ADAPTIVE_KEYWORDS, cs.RESTART, env.adap_ex_output,
            env.adap_ex_input, env.cpus, env.pele_exit_temp, env.residue, env.equil_steps, env.random_num)
        adaptive_exit.run(hook=True)
        env.logger.info("ExitPath Adaptive run successfully")


    if args.restart in ["all", "adaptive", "pele"]:

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
        center, radius, BS_sasa_min, BS_sasa_max = bx.create_box(args, env)
        env.logger.info("Box with center {} radius {} was created".format(center, radius))

        # Pele Exploration
        env.logger.info("Running standard Pele")
        ad.SimulationBuilder(env.pele_temp,  env.topology, cs.PELE_KEYWORDS, center, radius, BS_sasa_min, BS_sasa_max)
        adaptive_long = ad.SimulationBuilder(env.ad_l_temp,  env.topology, cs.ADAPTIVE_KEYWORDS,
            cs.RESTART, env.adap_l_output, env.adap_l_input, args.cpus, env.pele_temp, args.residue, env.random_num)
        adaptive_long.run()
        env.logger.info("Pele run successfully")

    if args.restart in ["all", "adaptive", "pele", "msm"]:

        # MSM Analysis
        env.logger.info("Running MSM analysis")
        msm.analyse_results(env, args, runTica=False)
        env.logger.info("MSM analysis run successfully")

        env.logger.info("{} System run successfully".format(args.residue))


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
    parser.add_argument("--test", action='store_true', help="Run a fast MSM_PELE test")
    parser.add_argument("--user_center", "-c", nargs='+', type=float, help='center of the box', default=None)
    parser.add_argument("--user_radius", "-r", type=float,  help="Radius of the box", default=None)
    parser.add_argument("--folder", "-wf", type=str,  help="Folder to apply the restart to", default=None)
    parser.add_argument("--pdb", action='store_true',  help="Use pdb files as output")
    
    args = parser.parse_args()
    if(args.clust > args.cpus and args.restart != "msm" and not args.test ):
        raise ValueError(cs.CLUSTER_ERROR.format(args.cpus, args.clust))
    else:
        run(args)
