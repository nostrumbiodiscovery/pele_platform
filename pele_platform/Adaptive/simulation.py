import os
import shutil
from AdaptivePELE import adaptiveSampling
import PPP.main as ppp

from pele_platform.Utilities.Helpers.map_atoms import AtomMapper
from pele_platform.Utilities.Helpers import helpers
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
from pele_platform.context import context


def run_adaptive():
    """
    Main function to prepare and launch the simulation.

    1) Prepare ligand and receptor
    2) Launch simulation
    3) Analyse simulation
    """
    missing_residues = []

    if context.parameters.adaptive_restart and not context.parameters.only_analysis:
        with helpers.cd(context.parameters.pele_dir):
            adaptiveSampling.main(context.parameters.ad_ex_temp)
            context.parameters.logger.info("Simulation run successfully (:\n\n")

    elif not context.parameters.only_analysis and not context.parameters.restart:
        context.parameters.logger.info(
            "System: {}; Platform Functionality: {}\n\n".format(
                context.parameters.residue, context.parameters.software))

        pdb_checker.PDBChecker(context.parameters.system).check_negative_residues()

        if context.parameters.perturbation:
            syst = sp.SystemBuilder.build_system(
                context.parameters.system,
                context.parameters.mae_lig,
                context.parameters.residue,
                context.parameters.pele_dir,
                inputs_dir=context.parameters.inputs_dir)
        else:
            syst = sp.SystemBuilder(
                context.parameters.system,
                None,
                None,
                context.parameters.pele_dir,
                inputs_dir=context.parameters.inputs_dir)

        context.parameters.logger.info("Prepare complex {}".format(syst.system))

        # If user overrides 'system' with 'input'
        if context.parameters.input:
            inputs_simulation = []

            for input_file in context.parameters.input:
                input_path = os.path.join(context.parameters.inputs_dir, os.path.basename(input_file))
                shutil.copy(input_file, input_path)

                if not context.parameters.no_ppp:
                    input_path, missing_residues, _, _, _ = ppp.main(
                        input_path,
                        context.parameters.inputs_dir,
                        output_pdb=["", ],
                        charge_terminals=context.parameters.charge_ter,
                        no_gaps_ter=context.parameters.gaps_ter,
                        constrain_smiles=None,
                        ligand_pdb=context.parameters.ligand_ref)

                input_path = os.path.join(context.parameters.inputs_dir, os.path.basename(input_path))
                inputs_simulation.append(input_path)

            context.parameters.input = inputs_simulation

        # If randomization in necessary (PPI, site_finder, global exploration)...

        elif context.parameters.full or context.parameters.randomize or context.parameters.ppi or context.parameters.site_finder:
            ligand_positions, box_radius, box_center = rd.randomize_starting_position(
                context.parameters.ligand_ref,
                context.parameters.receptor,
                outputfolder=context.parameters.inputs_dir,
                nposes=context.parameters.poses,
                test=context.parameters.test,
                user_center=context.parameters.center_of_interface,
                logger=context.parameters.logger)
            if not context.parameters.gpcr_orth and not context.parameters.out_in:
                context.parameters.box_center = box_center
                context.parameters.box_radius = box_radius

            if context.parameters.no_ppp:
                receptor = syst.system
            else:
                receptor, missing_residues, _, _, _ = ppp.main(
                    syst.system,
                    context.parameters.inputs_dir,  # to ensure it goes to pele_dir/input, not pele_dir
                    output_pdb=["", ],
                    charge_terminals=context.parameters.charge_ter,
                    no_gaps_ter=context.parameters.gaps_ter,
                    constrain_smiles=None,
                    ligand_pdb=context.parameters.ligand_ref)

            inputs = rd.join(
                receptor,
                ligand_positions,
                context.parameters.residue,
                output_folder=context.parameters.inputs_dir)

            context.parameters.input = [os.path.join(context.parameters.inputs_dir, inp)
                                for inp in inputs]
            context.parameters.system = context.parameters.input[0]

            hp.silentremove(ligand_positions)

        # Prepare System
        if context.parameters.no_ppp or context.parameters.input or context.parameters.covalent_residue:
            # No need to run system through PPP, if we already preprocessed
            # context.parameters.input
            if context.parameters.input:
                # If we have more than one input
                for input_file in context.parameters.input:
                    try:
                        shutil.copy(input_file, context.parameters.inputs_dir)
                    except shutil.SameFileError:  # systems that go through randomization are already moved
                        pass
            else:
                try:
                    shutil.copy(context.parameters.system, context.parameters.inputs_dir)
                except shutil.SameFileError:
                    pass
        else:
            context.parameters.nonstandard.extend(hp.find_nonstd_residue(syst.system))
            context.parameters.system, missing_residues, _, _, _ = ppp.main(
                syst.system,
                context.parameters.inputs_dir,
                output_pdb=["", ],
                charge_terminals=context.parameters.charge_ter,
                no_gaps_ter=context.parameters.gaps_ter,
                mid_chain_nonstd_residue=context.parameters.nonstandard,
                skip=context.parameters.skip_prep,
                back_constr=context.parameters.ca_constr,
                constrain_smiles=None,
                ligand_pdb=context.parameters.ligand_ref,
                ca_interval=context.parameters.ca_interval)

        context.parameters.constraints = alpha_constraints.retrieve_constraints(
            context.parameters.system,
            interval=context.parameters.ca_interval,
            back_constr=context.parameters.ca_constr,
            ter_constr=context.parameters.terminal_constr)

        # Metal constraints
        if not context.parameters.no_metal_constraints:
            metal_constraints, context.parameters.external_constraints = mc.main(
                context.parameters.external_constraints,
                context.parameters.input[0] if context.parameters.input else context.parameters.system,
                syst.system,
                permissive=context.parameters.permissive_metal_constr,
                all_metals=context.parameters.constrain_all_metals,
                external=context.parameters.external_constraints,
                logger=context.parameters.logger)

            context.parameters.external_constraints = hp.retrieve_constraints_for_pele(
                context.parameters.external_constraints, context.parameters.system)

            metal_constraints_json = hp.retrieve_constraints_for_pele(
                metal_constraints,
                context.parameters.system)

            context.parameters.external_constraints.extend(metal_constraints_json)

        else:
            context.parameters.external_constraints = hp.retrieve_constraints_for_pele(
                context.parameters.external_constraints, context.parameters.system)

        # Keep JSON ordered by having first title and then constraints
        if context.parameters.external_constraints:
            context.parameters.constraints = (context.parameters.constraints[0:1]
                                      + context.parameters.external_constraints
                                      + context.parameters.constraints[1:])
        if context.parameters.remove_constraints:
            context.parameters.constraints = ""

        context.parameters.logger.info(f"Complex {context.parameters.system} prepared\n\n")

        # Ligand/metal and solvent parameters
        if (context.parameters.perturbation or context.parameters.sidechain_perturbation) and context.parameters.use_peleffy:
            ligand_parametrizer = parametrizer.Parametrizer.from_parameters()
            ligand_parametrizer.parametrize_ligands_from(pdb_file=syst.system, ppp_file=context.parameters.system)

        elif (context.parameters.perturbation or context.parameters.sidechain_perturbation) and not context.parameters.use_peleffy:
            # Parametrize the ligand
            ligand_params = lg.LigandParametrization()
            ligand_params.generate()

            # Parametrize missing residues identified by PPP
            for res, __, _ in missing_residues:
                if res != context.parameters.residue and res not in context.parameters.skip_ligand_prep:
                    context.parameters.logger.info("Creating template for residue {}".format(res))
                    with hp.cd(context.parameters.pele_dir):
                        mr.create_template(res)
                    context.parameters.logger.info("Template {}z created\n\n".format(res))

        if context.parameters.ligand_conformations:
            ligand_conformations.LigandConformations(path=context.parameters.ligand_conformations, system=context.parameters.system,
                                                     resname=context.parameters.residue, forcefield=context.parameters.forcefield,
                                                     pele_dir=context.parameters.pele_dir).generate()

        # Create simulation box, if performing perturbation
        if context.parameters.perturbation:
            box = bx.BoxSetter(context.parameters.box_center,
                               context.parameters.box_radius,
                               context.parameters.ligand_ref,
                               context.parameters.logger)
            context.parameters.box = box.generate_json()
        else:
            context.parameters.box = ""

        # Solvent parameters
        solvent = sv.ImplicitSolvent(
            context.parameters.solvent,
            context.parameters.obc_tmp,
            context.parameters.template_folder,
            context.parameters.obc_file,
            context.parameters.logger,
        )
        solvent.generate()

        # Build PCA
        if context.parameters.pca_traj:
            pca_obj = pca.PCA(context.parameters.pca_traj, context.parameters.pele_dir)
            context.parameters.pca = pca_obj.generate(context.parameters.logger)

        # Core constraints based on SMILES string
        if context.parameters.constrain_core:
            smiles_input_pdb = context.parameters.input[0] if context.parameters.input else context.parameters.system
            smiles = smiles_constraints.SmilesConstraints(
                smiles_input_pdb,
                context.parameters.constrain_core,
                context.parameters.residue,
                context.parameters.chain,
                context.parameters.constrain_core_spring)
            smi_constraint = smiles.run()
            context.parameters.constraints = (context.parameters.constraints[0:1]
                                      + smi_constraint
                                      + context.parameters.constraints[1:])

        # Waters
        water_obj = wt.WaterIncluder(
            context.parameters.input if context.parameters.input else [context.parameters.system],
            context.parameters.n_waters,
            user_waters=context.parameters.waters,
            ligand_perturbation_params=context.parameters.parameters,
            water_center=context.parameters.water_center,
            water_radius=context.parameters.water_radius,
            allow_empty_selectors=context.parameters.allow_empty_selectors,
            water_temp=context.parameters.water_temp,
            water_trials=context.parameters.water_trials,
            water_overlap=context.parameters.water_overlap,
            water_constr=context.parameters.water_constr,
            test=context.parameters.test,
            water_freq=context.parameters.water_freq,
            ligand_residue=context.parameters.residue)
        water_obj.run()
        context.parameters.parameters = water_obj.ligand_perturbation_params
        context.parameters.water_ids_to_track = water_obj.water_ids_to_track

        AtomMapper(syst.system).run()

        # Metrics builder - builds JSON strings for PELE to be able to
        # track atom distances, RMSD, etc.
        metrics = mt.MetricBuilder()
        context.parameters.metrics = (
            metrics.distance_to_atom_json(
                os.path.join(
                    context.parameters.input[0] if context.parameters.input else context.parameters.system
                ),
                context.parameters.atom_dist,
            )
            if context.parameters.atom_dist
            else ""
        )
        context.parameters.native = (
            metrics.rsmd_to_json(context.parameters.native, context.parameters.chain) if context.parameters.native else ""
        )

        context.parameters.local_nonbonding_energy = metrics.local_nonbonding_energy_json(context.parameters.covalent_residue,
                                                                                  context.parameters.nonbonding_radius)

        # metal polarisation
        if context.parameters.polarize_metals:
            mp.change_metal_charges(
                context.parameters.template_folder,
                context.parameters.forcefield,
                context.parameters.polarization_factor,
                context.parameters.system,
            )

        context.parameters.adap_ex_input = ", ".join(
            [f'"{input_file}"' for input_file in context.parameters.input]).strip('"') if context.parameters.input else context.parameters.system

        # Fill in simulation templates
        adaptive = ad.SimulationBuilder(
            context.parameters.ad_ex_temp,
            context.parameters.pele_exit_temp,
            context.parameters.topology
        )
        adaptive.generate_inputs(water_obj)

        # Run simulation only if we are not in debug mode
        if not context.parameters.debug:
            context.parameters.logger.info("Running Simulation")
            adaptive.run()
            context.parameters.logger.info("Simulation run successfully (:\n\n")

    elif context.parameters.restart:
        # Start simulation from scratch (unlike adaptive_restart) but use files created in debug mode
        context.parameters.logger.info(f"Launching simulation from {context.parameters.pele_dir}")
        adaptive = ad.SimulationBuilder(context.parameters.ad_ex_temp, context.parameters.pele_exit_temp, context.parameters.topology)
        adaptive.run()
        context.parameters.logger.info("Simulation run successfully (:\n\n")

    # Run analysis
    if context.parameters.analyse and not context.parameters.debug:
        from pele_platform.analysis import Analysis

        # Retrieve water IDs to track from existing pele.conf, if running analysis only
        if context.parameters.only_analysis:
            context.parameters.water_ids_to_track = wt.water_ids_from_conf(context.parameters.pele_temp)

        analysis_folder = os.path.join(context.parameters.pele_dir, "results")

        analysis = Analysis.from_parameters()
        analysis.generate(
            analysis_folder,
            clustering_type=context.parameters.clustering_method.lower(),
            bandwidth=context.parameters.bandwidth,
            analysis_nclust=context.parameters.analysis_nclust,
            max_top_clusters=context.parameters.max_top_clusters,
            min_population=context.parameters.min_population,
            max_top_poses=context.parameters.max_top_poses,
            top_clusters_criterion=context.parameters.top_clusters_criterion,
            representatives_criterion=context.parameters.cluster_representatives_criterion)
