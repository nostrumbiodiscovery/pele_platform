import MSM_PELE.Utilities.Helpers.pele_env as pele
import MSM_PELE.constants as cs
import MSM_PELE.Utilities.Helpers.system_prep as sp
import MSM_PELE.Utilities.PPP.mut_prep4pele as ppp
import MSM_PELE.Utilities.PlopRotTemp.main as plop
import MSM_PELE.Utilities.Helpers.missing_residues as mr
import MSM_PELE.Utilities.Helpers.simulation as ad




def run_adaptive(args):
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
        ad.SimulationBuilder.simulation_handler(env, protein_constraints) 

    if args.restart in ["all", "glide_pele"] and args.software != "adaptive":
        # Run Adaptive Exit
        env.logger.info("Running Adaptive")
        adaptive_exit = ad.SimulationBuilder(env.ad_ex_temp, env.topology, cs.EX_ADAPTIVE_KEYWORDS, cs.RESTART, env.adap_ex_output,
            env.adap_ex_input, env.cpus, env.pele_exit_temp, env.residue, env.equil_steps, env.random_num)
        adaptive_exit.run()
        env.logger.info("Adaptive run successfully")

    elif args.software == "adaptive":
        env.logger.info("Running Adaptive")
        adaptive = ad.SimulationBuilder(args.adaptive, env.topology, cs.ADAPTIVE, env.adap_ex_input, env.cpus, args.pele, env.residue)
        ad.SimulationBuilder(args.pele,  env.topology, ["CHAIN"], env.chain)
        adaptive.run()
        env.logger.info("Adaptive run successfully")




    return env
