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
import MSM_PELE.Helpers.constraints as ct
import MSM_PELE.Helpers.system_prep as sp
import MSM_PELE.Helpers.box as bx
import MSM_PELE.PPP.mut_prep4pele as ppp
import MSM_PELE.Helpers.msm_analysis as msm
import MSM_PELE.Helpers.missing_residues as mr


def run(system, residue, chain, mae_lig, user_box, charge_ter, gaps_ter, clusters, forcefield, confile, native, cpus, core, mtor, n, clean, restart, gridres):

    # Build folders and logging
    env = pele.Pele_env_Builder(cs.FOLDERS, cs.FILES, forcefield, system, residue, cpus, restart, native, chain, mae_lig)
    env.create()

    if restart == "all":

        # Building system and ligand
        env.logger.info("Checking {} system for Pele".format(residue))
        if mae_lig:
            receptor = system
            lig = mae_lig
            lig_ref = sp.convert_pdb(mae_lig, env.pele_dir)
            system = sp.build_complex(receptor, lig_ref, env.pele_dir)
        else:
            receptor, lig_ref = sp.retrieve_receptor(system, residue, env.pele_dir)
            lig, residue = sp.convert_mae(lig_ref)
        system_fix, missing_residues, gaps, metals = ppp.main(system, env.pele_dir, charge_terminals=charge_ter, no_gaps_ter=gaps_ter)
        protein_constraints = ct.retrieve_constraints(system_fix, gaps, metals)
        env.logger.info(cs.SYSTEM.format(system_fix, missing_residues, gaps, metals))

        # Parametrize Ligands
        env.logger.info("Creating template for residue {}".format(residue))
        if mae_lig:
            mae_charges = True
            template, rotamers_file = plop.main(mae_lig, residue, env.pele_dir, forcefield, mtor, n, core, mae_charges, clean, gridres)
            hp.silentremove([system])
        else:
            mae_charges = False
            template, rotamers_file = plop.main(lig, residue, env.pele_dir, forcefield, mtor, n, core, mae_charges, clean, gridres)
            hp.silentremove([lig])
        env.logger.info("Template {} created".format(template))

        for res, __, _ in missing_residues:
            if res != residue:
                env.logger.info("Creating template for residue {}".format(res))
                mr.create_template(system_fix, res, env.pele_dir, forcefield)
                env.logger.info("Template {}z created".format(res))

        # Fill in Templates
        ad.SimulationBuilder(env.pele_exit_temp, cs.EX_PELE_KEYWORDS, native, forcefield, chain, "\n".join(protein_constraints), cpus, env.license)
        ad.SimulationBuilder(env.pele_temp, cs.EX_PELE_KEYWORDS, native, forcefield, chain, "\n".join(protein_constraints), cpus, env.license)

    if restart in ["all", "adaptive"]:
        # Adaptive Exit
        env.logger.info("Running ExitPath Adaptive")
        adaptive_exit = ad.SimulationBuilder(env.ad_ex_temp, cs.EX_ADAPTIVE_KEYWORDS, cs.RESTART, env.adap_ex_output,
            env.adap_ex_input, env.cpus, env.pele_exit_temp, env.residue, env.equil_steps, env.random_num)
        adaptive_exit.run(hook=True)
        env.logger.info("ExitPath Adaptive run successfully")


    if restart in ["all", "pele"]:

        # KMeans Clustering
        if not os.path.isfile(env.clusters_output):
            env.logger.info("Running MSM Clustering")
            with hp.cd(env.adap_ex_output):
                cl.main(num_clusters=clusters, output_folder=env.cluster_output, ligand_resname=residue, atom_ids="")
            env.logger.info("MSM Clustering run successfully")
        else:
            pass

        # Create Box
        env.logger.info("Creating box")
        bx.is_exit_finish(env.adap_ex_output)
        if not user_box:
            center_mass = cm.center_of_mass(env.ligand_ref)
            center, radius = bx.main(env.adap_ex_input, env.clusters_output, center_mass)
            bx.build_box(center, radius, env.box_temp)
        else:
            center, radius = bx.retrieve_box_info(user_box, env.clusters_output)
            shutil.copy(user_box, os.path.join(env.pele_dir, "box.pdb"))
        env.logger.info("Box with center {} radius {} was created".format(center, radius))

        # Pele Exploration
        env.logger.info("Running standard Pele")
        ad.SimulationBuilder(env.pele_temp, cs.PELE_KEYWORDS, center, radius)
        adaptive_long = ad.SimulationBuilder(env.ad_l_temp, cs.ADAPTIVE_KEYWORDS,
            cs.RESTART, env.adap_l_output, env.adap_l_input, cpus, env.pele_temp, residue, env.random_num)
        adaptive_long.run()
        env.logger.info("Pele run successfully")

    if restart in ["all", "pele", "msm"]:

        # MSM Analysis
        env.logger.info("Running MSM analysis")
        msm.analyse_results(env.adap_l_output, residue, cpus, env.pele_dir)
        env.logger.info("MSM analysis run successfully")

        env.logger.info("{} System run successfully".format(residue))


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Run Adaptive Pele Platform')
    parser.add_argument('input', type=str, help='complex to run pele on')
    parser.add_argument('residue', type=str, help='residue of the ligand to extract', default=cs.LIG_RES)
    parser.add_argument('chain', type=str, help='chain of the ligand to extract', default=cs.LIG_CHAIN)
    parser.add_argument("--mae_lig", type=str, help="ligand .mae file to include QM charges coming from jaguar")
    parser.add_argument("--box", type=str, help="Exploration box for Pele")
    parser.add_argument("--charge_ter", help="Charge protein terminals", action='store_true')
    parser.add_argument("--gaps_ter", help="Include TER when a possible gap is found", action='store_true')
    parser.add_argument("--clust", type=int, help="Numbers of clusters to start PELE's exploration with", default=cs.CLUSTERS)
    parser.add_argument('--forc', type=str, help='chain of the ligand to extract', default=cs.FORCEFIELD)
    parser.add_argument('--confile', type=str, help='your own pele configuration file', default=cs.PELE_CONFILE)
    parser.add_argument('--native', type=str, help='native file to compare RMSD to', default="")
    parser.add_argument('--cpus', type=int, help='number of processors', default=cs.CPUS)
    parser.add_argument("--core", type=int, help="Give one atom of the core section", default=-1)
    parser.add_argument("--mtor", type=int, help="Gives the maximum number of torsions allowed in each group.  Will freeze bonds to extend the core if necessary.", default=4)
    parser.add_argument("--n", type=int, help="Maximum Number of Entries in Rotamer File", default=1000)
    parser.add_argument("--clean", help="Whether to clean up all the intermediate files", action='store_true')
    parser.add_argument("--restart", type=str, help="Restart the platform from [all, pele, msm] with these keywords", default=cs.PLATFORM_RESTART)
    parser.add_argument("--gridres", type=str, help="Rotamers angle resolution", default=cs.GRIDRES)
    args = parser.parse_args()
    print(args.restart)
    if(args.clust > args.cpus and args.restart != "msm"):
        raise ValueError(cs.CLUSTER_ERROR.format(args.cpus, args.clust))
    else:
        run(args.input, args.residue, args.chain, args.mae_lig, args.box, args.charge_ter, args.gaps_ter, args.clust, args.forc, args.confile, args.native, args.cpus, args.core, args.mtor, args.n, args.clean, args.restart, args.gridres)
