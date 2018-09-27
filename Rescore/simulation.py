import subprocess
import PELE_Platform.Utilities.Helpers.pele_env as pele
import PELE_Platform.constants as cs
import PELE_Platform.Utilities.Helpers.system_prep as sp
import PELE_Platform.Utilities.PPP.mut_prep4pele as ppp
import PELE_Platform.Utilities.PlopRotTemp.main as plop
import PELE_Platform.Utilities.Helpers.missing_residues as mr
import PELE_Platform.Utilities.Helpers.simulation as ad
import PELE_Platform.Utilities.Helpers.center_of_mass as cm
import PELE_Platform.Utilities.Helpers.randomize as rd
import PELE_Platform.Utilities.Helpers.helpers as hp



def run_adaptive(args):
    # Build folders and logging
    env = pele.EnviroBuilder.build_env(args)

    if args.restart == "all":

        # Build System
        env.logger.info("Checking {} system for Pele".format(args.residue))
        syst = sp.SystemBuilder.build_system(args.system, args.mae_lig, args.residue, env.pele_dir)
        if args.full:
            ligand_positions = rd.randomize_starting_position(env.ligand_ref, "input_ligand.pdb", env.residue, env.receptor, None, None, env)
            inputs = rd.join(env.receptor, ligand_positions, env)
            hp.silentremove(ligand_positions)
            #Parsing input for errors and saving them as inputs
            env.adap_ex_input = ", ".join([ '"' + ppp.main(input, env.pele_dir, output_pdb=["" , ], 
                charge_terminals=args.charge_ter, no_gaps_ter=args.gaps_ter)[0] + '"' for input in inputs ])
            hp.silentremove(inputs)
        # Prepare System
        system_fix, missing_residues, gaps, metals, env.constraints = ppp.main(syst.system, env.pele_dir, charge_terminals=args.charge_ter, no_gaps_ter=args.gaps_ter)
        env.logger.info(cs.SYSTEM.format(system_fix, missing_residues, gaps, metals))

        # Parametrize Ligand
        if not env.external_template and not env.external_rotamers:
            env.logger.info("Creating template for residue {}".format(args.residue))
            plop.parametrize_miss_residues(args, env, syst)
            env.logger.info("Template {}z created".format(args.residue.lower()))
        else:
            cmd_to_move_template = "cp {} {}".format(env.external_template,  env.template_folder)
            subprocess.call(cmd_to_move_template.split())
            cmd_to_move_rotamer_file = "cp {} {}".format(env.external_rotamers,  env.rotamers_folder)
            subprocess.call(cmd_to_move_rotamer_file.split())



        # Parametrize missing residues
        for res, __, _ in missing_residues:
            if res != args.residue:
                env.logger.info("Creating template for residue {}".format(res))
                mr.create_template(system_fix, res, env.pele_dir, args.forcefield)
                env.logger.info("Template {}z created".format(res))

        #Set Box
        env.box_center = cm.center_of_mass(env.ligand_ref)

        # Fill in Simulation Templates
        adaptive = ad.SimulationBuilder(env.ad_ex_temp, env.pele_exit_temp, env)
        adaptive.run()
        
    return env
