import os
import shutil
from AdaptivePELE import adaptiveSampling
import PPP.main as ppp

from pele_platform.Utilities.Helpers.map_atoms import AtomMapper
from pele_platform.Utilities.Helpers import helpers
from pele_platform.Utilities.Parameters.parameters import ParametersBuilder
from pele_platform.Utilities.Helpers.constraints import (
    alpha_constraints,
    smiles_constraints,
)
import pele_platform.Utilities.Helpers.simulation as ad
import pele_platform.Utilities.Helpers.system_prep as sp
import pele_platform.Utilities.Helpers.missing_residues as mr
import pele_platform.Utilities.Helpers.randomize as rd
import pele_platform.Utilities.Helpers.helpers as hp
import pele_platform.Utilities.Helpers.Metals.metal_constraints as mc
import pele_platform.Utilities.Helpers.Metals.metal_polarisation as mp
import pele_platform.Adaptive.metrics as mt
import pele_platform.Utilities.Helpers.water as wt
from pele_platform.Adaptive import parametrizer
import pele_platform.Adaptive.box as bx
import pele_platform.Adaptive.pca as pca
import pele_platform.Adaptive.plop_solvent as sv
import pele_platform.Adaptive.plop_ligand_parametrization as lg
from pele_platform.Utilities.Helpers import ligand_conformations
from pele_platform.Checker import pdb_checker


def run_adaptive(args):
    """
    Main function to prepare and launch the simulation.

    1) Create EnviroBuilder - variables, folders, logger...
    2) Prepare ligand and receptor
    3) Launch simulation
    4) Analyse simulation
    """
    builder = ParametersBuilder()
    parameters = builder.build_adaptive_variables(args)
    parameters.create_files_and_folders()
    shutil.copy(args.yamlfile, parameters.pele_dir)
    missing_residues = []

    if parameters.adaptive_restart and not parameters.only_analysis:
        with helpers.cd(parameters.pele_dir):
            adaptiveSampling.main(parameters.ad_ex_temp)
            parameters.logger.info("Simulation run successfully (:\n\n")

    elif not parameters.only_analysis and not parameters.restart:
        parameters.logger.info(
            "System: {}; Platform Functionality: {}\n\n".format(
                parameters.residue, parameters.software))

        pdb_checker.PDBChecker(parameters.system).check_negative_residues()

        # Create inputs directory
        if not os.path.exists(parameters.inputs_dir):
            os.mkdir(parameters.inputs_dir)

        if parameters.perturbation:
            syst = sp.SystemBuilder.build_system(
                parameters.system,
                args.mae_lig,
                args.residue,
                parameters.pele_dir,
                inputs_dir=parameters.inputs_dir)
        else:
            syst = sp.SystemBuilder(
                parameters.system,
                None,
                None,
                parameters.pele_dir,
                inputs_dir=parameters.inputs_dir)

        parameters.logger.info("Prepare complex {}".format(syst.system))

        # If user overrides 'system' with 'input'
        if parameters.input:
            parameters.inputs_simulation = []

            for input in parameters.input:
                input_path = os.path.join(
                    parameters.inputs_dir, os.path.basename(input))
                shutil.copy(input, input_path)

                if parameters.no_ppp:
                    input_proc = input_path
                else:
                    input_proc, missing_residues, _, _, _ = ppp.main(
                            input_path,
                            parameters.inputs_dir,  # to ensure it goes to pele_dir/inputs, not pele_dir
                            output_pdb=["", ],
                            charge_terminals=args.charge_ter,
                            no_gaps_ter=args.gaps_ter,
                            constrain_smiles=None,
                            ligand_pdb=parameters.ligand_ref)
                input_proc = os.path.basename(input_proc)
                input_proc = os.path.join(parameters.inputs_dir, input_proc)
                parameters.inputs_simulation.append(input_proc)
            parameters.adap_ex_input = ", ".join(
                ['"' + input + '"' for input in parameters.inputs_simulation]
            ).strip('"')

        # If randomization in necessary (PPI, site_finder, global exploration)...
        elif args.full or args.randomize or args.ppi or args.site_finder:
            ligand_positions, box_radius, box_center = rd.randomize_starting_position(
                parameters.ligand_ref,
                parameters.receptor,
                outputfolder=parameters.inputs_dir,
                nposes=parameters.poses,
                test=parameters.test,
                user_center=parameters.center_of_interface,
                logger=parameters.logger)
            if not args.gpcr_orth and not args.out_in:
                parameters.box_center = box_center
                parameters.box_radius = box_radius

            if parameters.no_ppp:
                receptor = syst.system
            else:
                receptor, missing_residues, _, _, _ = ppp.main(
                    syst.system,
                    parameters.inputs_dir,  # to ensure it goes to pele_dir/input, not pele_dir
                    output_pdb=["", ],
                    charge_terminals=args.charge_ter,
                    no_gaps_ter=args.gaps_ter,
                    constrain_smiles=None,
                    ligand_pdb=parameters.ligand_ref)

            inputs = rd.join(
                receptor,
                ligand_positions,
                parameters.residue,
                output_folder=parameters.inputs_dir)

            parameters.input = [os.path.join(parameters.inputs_dir, inp)
                                for inp in inputs]

            parameters.adap_ex_input = ", ".join(
                ['"' + input + '"' for input in parameters.input]
            ).strip('"')
            hp.silentremove(ligand_positions)

        # Prepare System
        if parameters.no_ppp or parameters.input:  # No need to run system through PPP, if we already preprocessed
            # parameters.input
            if parameters.input:
                # If we have more than one input
                for input in parameters.input:
                    try:
                        shutil.copy(input, parameters.inputs_dir)
                    except shutil.SameFileError:  # systems that go through randomization are already moved
                        pass
            else:
                shutil.copy(parameters.system, parameters.inputs_dir)
        else:
            parameters.nonstandard.extend(hp.find_nonstd_residue(syst.system))
            parameters.system, missing_residues, _, _, _ = ppp.main(
                syst.system,
                parameters.inputs_dir,
                output_pdb=["", ],
                charge_terminals=args.charge_ter,
                no_gaps_ter=args.gaps_ter,
                mid_chain_nonstd_residue=parameters.nonstandard,
                skip=parameters.skip_prep,
                back_constr=parameters.ca_constr,
                constrain_smiles=None,
                ligand_pdb=parameters.ligand_ref,
                ca_interval=parameters.ca_interval)

        parameters.constraints = alpha_constraints.retrieve_constraints(
            parameters.system,
            interval=parameters.ca_interval,
            back_constr=parameters.ca_constr,
            ter_constr=parameters.terminal_constr)

        # Metal constraints
        if not args.no_metal_constraints:
            metal_constraints, parameters.external_constraints = mc.main(
                args.external_constraints,
                os.path.join(
                    parameters.inputs_dir,
                    parameters.adap_ex_input.split(",")[0].strip().strip('"')),
                syst.system,
                permissive=parameters.permissive_metal_constr,
                all_metals=args.constrain_all_metals,
                external=parameters.external_constraints,
                logger=parameters.logger)

            parameters.external_constraints = hp.retrieve_constraints_for_pele(
                parameters.external_constraints, parameters.system)

            metal_constraints_json = hp.retrieve_constraints_for_pele(
                metal_constraints,
                os.path.join(
                    parameters.inputs_dir,
                    parameters.adap_ex_input.split(",")[0].strip().strip('"')))

            parameters.external_constraints.extend(metal_constraints_json)

        else:
            parameters.external_constraints = hp.retrieve_constraints_for_pele(
                parameters.external_constraints, parameters.system)

        # Keep JSON ordered by having first title and then constraints
        if parameters.external_constraints:
            parameters.constraints = (parameters.constraints[0:1]
                                      + parameters.external_constraints
                                      + parameters.constraints[1:])
        if parameters.remove_constraints:
            parameters.constraints = ""

        parameters.logger.info(f"Complex {parameters.system} prepared\n\n")

        # Ligand/metal and solvent parameters
        if (parameters.perturbation or parameters.sidechain_perturbation) and parameters.use_peleffy:
            ligand_parametrizer = parametrizer.Parametrizer.from_parameters(parameters)
            ligand_parametrizer.parametrize_ligands_from(pdb_file=syst.system, ppp_file=parameters.system)

        elif (parameters.perturbation or parameters.sidechain_perturbation) and not parameters.use_peleffy:
            # Parametrize the ligand
            ligand_params = lg.LigandParametrization(parameters)
            ligand_params.generate()

            # Parametrize missing residues identified by PPP
            for res, __, _ in missing_residues:
                if res != args.residue and res not in parameters.skip_ligand_prep:
                    parameters.logger.info("Creating template for residue {}".format(res))
                    with hp.cd(parameters.pele_dir):
                        mr.create_template(parameters, res)
                    parameters.logger.info("Template {}z created\n\n".format(res))

        # Covalent residue parametrization should not run in refinement simulation
        if parameters.covalent_residue and os.path.basename(parameters.pele_dir) != "2_refinement":
            parametrizer.parametrize_covalent_residue(parameters.pele_data, parameters.pele_dir, parameters.gridres,
                                                      parameters.residue_type, parameters.residue,
                                                      ppp_system=parameters.system)

        if parameters.ligand_conformations:
            ligand_conformations.LigandConformations(path=parameters.ligand_conformations, system=parameters.system,
                                                     resname=parameters.residue, forcefield=parameters.forcefield,
                                                     pele_dir=parameters.pele_dir).generate()

        # Create simulation box, if performing perturbation
        if parameters.perturbation:
            box = bx.BoxSetter(parameters.box_center,
                               parameters.box_radius,
                               parameters.ligand_ref,
                               parameters.logger)
            parameters.box = box.generate_json()
        else:
            parameters.box = ""

        # Solvent parameters
        solvent = sv.ImplicitSolvent(
            parameters.solvent,
            parameters.obc_tmp,
            parameters.template_folder,
            parameters.obc_file,
            parameters.logger,
        )
        solvent.generate()

        # Build PCA
        if parameters.pca_traj:
            pca_obj = pca.PCA(parameters.pca_traj, parameters.pele_dir)
            parameters.pca = pca_obj.generate(parameters.logger)

        # Core constraints based on SMILES string
        if parameters.constrain_core:
            smiles_input_pdb = os.path.join(
                parameters.inputs_dir, parameters.adap_ex_input.split(",")[0]
            )
            smiles = smiles_constraints.SmilesConstraints(
                smiles_input_pdb,
                parameters.constrain_core,
                parameters.residue,
                parameters.chain,
                parameters.constrain_core_spring)
            smi_constraint = smiles.run()
            parameters.constraints = (parameters.constraints[0:1]
                                      + smi_constraint
                                      + parameters.constraints[1:])

        # Waters
        input_waters = [input.strip().strip('"')
                        for input in parameters.adap_ex_input.split(",")]
        input_waters = [os.path.join(parameters.inputs_dir, inp)
                        for inp in input_waters]
        water_obj = wt.WaterIncluder(
            input_waters,
            parameters.n_waters,
            user_waters=parameters.waters,
            ligand_perturbation_params=parameters.parameters,
            water_center=args.water_center,
            water_radius=parameters.water_radius,
            allow_empty_selectors=parameters.allow_empty_selectors,
            water_temp=parameters.water_temp,
            water_trials=parameters.water_trials,
            water_overlap=parameters.water_overlap,
            water_constr=parameters.water_constr,
            test=parameters.test,
            water_freq=parameters.water_freq,
            ligand_residue=parameters.residue)
        water_obj.run()
        parameters.parameters = water_obj.ligand_perturbation_params
        parameters.water_ids_to_track = water_obj.water_ids_to_track

        # Check if atoms need mapping due to preprocessing
        args = AtomMapper(args, parameters, syst.system).run()

        # Metrics builder - builds JSON strings for PELE to be able to
        # track atom distances, RMSD, etc.
        metrics = mt.MetricBuilder()
        parameters.metrics = (
            metrics.distance_to_atom_json(
                os.path.join(
                    parameters.inputs_dir,
                    parameters.adap_ex_input.split(",")[0].strip().strip('"'),
                ),
                args.atom_dist,
            )
            if args.atom_dist
            else ""
        )
        parameters.native = (
            metrics.rmsd_to_json(args.native, parameters.chain) if args.native else ""
        )

        parameters.local_nonbonding_energy = metrics.local_nonbonding_energy_json(parameters.covalent_residue,
                                                                                  parameters.nonbonding_radius)

        # metal polarisation
        if parameters.polarize_metals:
            mp.change_metal_charges(
                parameters.template_folder,
                parameters.forcefield,
                parameters.polarization_factor,
                parameters.system,
            )

        # Point adaptive.conf to input dir
        parameters.adap_ex_input = os.path.join(
            parameters.inputs_dir, parameters.adap_ex_input
        )

        # Fill in simulation templates
        adaptive = ad.SimulationBuilder(
            parameters.ad_ex_temp, parameters.pele_exit_temp, parameters.topology
        )
        adaptive.generate_inputs(parameters, water_obj)

        # Run simulation only if we are not in debug mode
        if not parameters.debug:
            parameters.logger.info("Running Simulation")
            adaptive.run()
            parameters.logger.info("Simulation run successfully (:\n\n")

    elif parameters.restart:
        # Start simulation from scratch (unlike adaptive_restart) but use files created in debug mode
        parameters.logger.info(f"Launching simulation from {parameters.pele_dir}")
        adaptive = ad.SimulationBuilder(parameters.ad_ex_temp, parameters.pele_exit_temp, parameters.topology)
        adaptive.run()
        parameters.logger.info("Simulation run successfully (:\n\n")

    # Run analysis
    if parameters.analyse and not parameters.debug:
        from pele_platform.analysis import Analysis

        # Retrieve water IDs to track from existing pele.conf, if running analysis only
        if parameters.only_analysis:
            parameters.water_ids_to_track = wt.water_ids_from_conf(parameters.pele_temp)

        analysis_folder = os.path.join(parameters.pele_dir, "results")

        analysis = Analysis.from_parameters(parameters)
        analysis.generate(
            analysis_folder,
            clustering_type=parameters.clustering_method.lower(),
            bandwidth=parameters.bandwidth,
            analysis_nclust=parameters.analysis_nclust,
            max_top_clusters=parameters.max_top_clusters,
            min_population=parameters.min_population,
            max_top_poses=parameters.max_top_poses,
            top_clusters_criterion=parameters.top_clusters_criterion,
            representatives_criterion=parameters.cluster_representatives_criterion)

    return parameters
