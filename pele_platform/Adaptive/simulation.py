import os
import shutil
import AdaptivePELE.adaptiveSampling as adt
import PPP.main as ppp
from pele_platform.Utilities.Helpers import helpers
import pele_platform.Utilities.Parameters.pele_env as pele
import pele_platform.Utilities.Helpers.constraints as ct
import pele_platform.constants.constants as cs
import pele_platform.Utilities.Helpers.simulation as ad
import pele_platform.Utilities.Helpers.system_prep as sp
import pele_platform.Utilities.Helpers.missing_residues as mr
import pele_platform.Utilities.Helpers.randomize as rd
import pele_platform.Utilities.Helpers.helpers as hp
import pele_platform.Utilities.Helpers.Metals.metal_constraints as mc
import pele_platform.Utilities.Helpers.Metals.metal_polarisation as mp
import pele_platform.Adaptive.metrics as mt
import pele_platform.Utilities.Helpers.water as wt
import pele_platform.Analysis.plots as pt
import pele_platform.Utilities.Parameters.pele_env as pv
import pele_platform.Adaptive.ligand_parametrization as lg
import pele_platform.Adaptive.box as bx
import pele_platform.Adaptive.solvent as sv
import pele_platform.Adaptive.pca as pca
# import pele_platform.RNA.prep as pr


def run_adaptive(args: pv.EnviroBuilder) -> pv.EnviroBuilder:
    """
    Main function to prepare and launch simulation

    1) Crate working folders and inputs
    2) Prepare ligand and receptor
    3) Launch simulation
    4) Analyse simulation
    """

    # Build Folders and Logging and env variable that will contain
    # all main attributes of the simulation
    env = pele.EnviroBuilder()
    env.software = "Adaptive"
    env.build_adaptive_variables(args)
    env.create_files_and_folders()
    shutil.copy(args.yamlfile, env.pele_dir)

    if env.adaptive_restart and not env.only_analysis:
        with helpers.cd(env.pele_dir):
            adt.main(env.ad_ex_temp)
            env.logger.info("Simulation run successfully (:\n\n")

    elif not env.only_analysis:

        env.logger.info("System: {}; Platform Functionality: {}\n\n".format(env.residue, env.software))
        
        if env.perturbation:
            syst = sp.SystemBuilder.build_system(env.system, args.mae_lig, args.residue, env.pele_dir)
        else:
            syst = sp.SystemBuilder(env.system, None, None, env.pele_dir)
        
        env.logger.info("Prepare complex {}".format(syst.system))
           
        # Choose your own input
        
        # User specifies more than one input
        if env.input:
            env.inputs_simulation = []
            for input in env.input:
                input_path  = os.path.join(env.pele_dir, os.path.basename(input))
                shutil.copy(input, input_path)
                
                if env.no_ppp:
                    input_proc = input
                else:
                    input_proc = os.path.basename(ppp.main(input_path, env.pele_dir, output_pdb=["" , ],
                                charge_terminals=args.charge_ter, no_gaps_ter=args.gaps_ter,
                                constrain_smiles=env.constrain_smiles, ligand_pdb=env.ligand_ref)[0])
                env.inputs_simulation.append(input_proc)
                hp.silentremove([input_path])
            env.adap_ex_input = ", ".join(['"' + input + '"' for input in env.inputs_simulation]).strip('"')
        
        # Global exploration mode: Create inputs around protein
        elif args.full or args.randomize or args.ppi:
            ligand_positions, box_radius, box_center = rd.randomize_starting_position(env.ligand_ref, env.receptor,
                outputfolder=env.pele_dir, nposes=env.poses, test=env.test, user_center=env.center_of_interface, logger=env.logger)
            if not args.gpcr_orth:
                env.box_center = box_center
                env.box_radius = box_radius
            if env.no_ppp:
                receptor = syst.system
            else:
                receptor = ppp.main(syst.system, env.pele_dir, output_pdb=["" , ],
                            charge_terminals=args.charge_ter, no_gaps_ter=args.gaps_ter,
                            constrain_smiles=env.constrain_smiles, ligand_pdb=env.ligand_ref)[0]
            inputs = rd.join(receptor, ligand_positions, env.residue, output_folder=env.pele_dir)
            env.adap_ex_input = ", ".join(['"' + os.path.basename(input) + '"' for input in inputs]).strip('"')
            hp.silentremove(ligand_positions)
            # Parsing input for errors and saving them as inputs

        # Prepare System
        if env.no_ppp:
            missing_residues = []
            gaps = {}
            metals = {}
            env.constraints = ct.retrieve_constraints(env.system, gaps, metals, back_constr=env.ca_constr)
            if env.input:
                # If we have more than one input
                for input in env.input:
                    shutil.copy(input, env.pele_dir)
            else:
                shutil.copy(env.system, env.pele_dir)
        else:
            env.nonstandard.extend(hp.find_nonstd_residue(syst.system))
            env.system, missing_residues, gaps, metals, env.constraints = ppp.main(syst.system, env.pele_dir, output_pdb=["" , ], charge_terminals=args.charge_ter, no_gaps_ter=args.gaps_ter, mid_chain_nonstd_residue=env.nonstandard, skip=env.skip_prep, back_constr=env.ca_constr, constrain_smiles=env.constrain_smiles, ligand_pdb=env.ligand_ref)

        # Metal constraints
        if not args.no_metal_constraints:
            metal_constraints, env.external_constraints = mc.main(args.external_constraints, os.path.join(env.pele_dir, env.adap_ex_input.split(",")[0].strip().strip('"')), syst.system, permissive=env.permissive_metal_constr, all_metals=args.constrain_all_metals, external=env.external_constraints, logger=env.logger)
            env.external_constraints = hp.retrieve_constraints_for_pele(env.external_constraints, env.system)
            metal_constraints_json = hp.retrieve_constraints_for_pele(metal_constraints, env.system)
            env.external_constraints.extend(metal_constraints_json)
        else:
            env.external_constraints = hp.retrieve_constraints_for_pele(env.external_constraints, env.system)

        # Keep JSON ordered by having first title and then constraints
        if env.external_constraints:
            env.constraints = env.constraints[0:1] + env.external_constraints + env.constraints[1:]
        if env.remove_constraints:
            env.constraints = ""
        env.logger.info(cs.SYSTEM.format(missing_residues, gaps, metals))
        env.logger.info("Complex {} prepared\n\n".format(env.system))

        # Ligand parameters and simulation box
        if env.perturbation:
            ligand_params = lg.LigandParametrization(env)
            ligand_params.generate()
            box = bx.BoxSetter(env.box_center, env.box_radius, env.ligand_ref, env.logger)
            env.box = box.generate_json()
        else:
            env.box = ""

        # Parametrize missing residues
        for res, __, _ in missing_residues:
            if res != args.residue and res not in env.skip_ligand_prep:
                env.logger.info("Creating template for residue {}".format(res))
                with hp.cd(env.pele_dir):
                    mr.create_template(env, res)
                env.logger.info("Template {}z created\n\n".format(res))
    
        # Solvent parameters
        solvent = sv.ImplicitSolvent(env.solvent, env.obc_tmp, env.template_folder, env.obc_file, env.logger)
        solvent.generate()
        
        # Build PCA
        if env.pca_traj:
            pca_obj = pca.PCA(env.pca_traj, env.pele_dir)
            env.pca = pca_obj.generate(env.logger)
        
        ####### Add waters, if needed
        input_waters = [input.strip().strip('"') for input in env.adap_ex_input.split(",")]
        input_waters = [os.path.join(env.pele_dir, inp) for inp in input_waters]
        water_obj = wt.WaterIncluder(input_waters, env.n_waters, 
            user_waters=env.waters, ligand_perturbation_params=env.parameters, 
            water_center=args.water_center, water_radius=env.water_radius,
            allow_empty_selectors=env.allow_empty_selectors, water_temp=env.water_temp,
            water_trials=env.water_trials, water_overlap=env.water_overlap,
            water_constr=env.water_constr, test=env.test, water_freq=env.water_freq,
            ligand_residue=env.residue)
        water_obj.run()
        env.parameters = water_obj.ligand_perturbation_params

        #metrics
        metrics = mt.MetricBuilder()
        if args.atom_dist:
            env.metrics = metrics.distance_to_atom_json(env.system, args.atom_dist)
        else:
            env.metrics = ""
        if args.native:
            env.native = metrics.rsmd_to_json(args.native, env.chain)
        else:
            env.native = ""

        # metal polarisation
        if env.polarize_metals:
            mp.change_metal_charges(env.template_folder, env.forcefield, env.polarization_factor, env.system)

        ############Fill in Simulation Templates############
        # Fill in simulation templates
        adaptive = ad.SimulationBuilder(env.ad_ex_temp,
            env.pele_exit_temp, env.topology)
        adaptive.generate_inputs(env, water_obj)

    if env.analyse and not env.debug:
        report = pt.analyse_simulation(env.report_name, env.traj_name[:-4]+"_", 
            os.path.join(env.pele_dir, env.output), env.residue, cpus=env.cpus,
            output_folder=env.pele_dir, clustering=env.perturbation, mae=env.mae,
            nclusts=env.analysis_nclust, overwrite=env.overwrite, topology=env.topology,
            be_column=env.be_column, limit_column=env.limit_column, te_column=env.te_column, logger=env.logger)
        env.logger.info("Pdf summary report successfully written to: {}".format(report))
    return env
