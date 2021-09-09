from abc import abstractmethod
import glob
import os
import re

import pele_platform.Adaptive.simulation as si
from pele_platform.Errors import custom_errors
from pele_platform.Utilities.Helpers import helpers
from pele_platform.building_blocks.preparation import prepare_structure
from pele_platform.building_blocks import blocks
from pele_platform.Utilities.Parameters import parameters
import pele_platform.features.adaptive as ft
from pele_platform.Adaptive import parametrizer
from pele_platform.context import context


class Simulation(blocks.Block):
    """
    Base Simulation class to run all simulation types.

    One class to rule them all, one class to find them, one class to bring them all and in PELE bind them.
    """
    keyword = None

    def run(self) -> (parameters.ParametersBuilder, parameters.Parameters):
        self.simulation_setup()
        si.run_adaptive()
        self.post_simulation_cleanup()

    def simulation_setup(self):
        self.set_working_folder()
        self.set_user_params()
        self.restart_checker()
        self.create_folders()
        self.set_params(simulation_type=self.keyword)
        self.set_package_params()
        if hasattr(context.parameters, "next_step"):
            context.parameters.input = glob.glob(context.parameters.next_step)
        self.water_handler()

    def post_simulation_cleanup(self):
        """
        Sets both types of restart to false, so that the platform does not look for logger file while executing the
        next building block. Disables PDB preprocessing, since it would have been done at the first step.
        """
        context.parameters.restart = False
        context.parameters.adaptive_restart = False
        context.parameters.no_ppp = True

    @abstractmethod
    def set_package_params(self):
        pass

    def restart_checker(self):
        """
        Check what kind of restart should be executed.
        If the output folder exists and the user wants to restart adaptive -> adaptive_restart.
        If the output folder does not exists and the user wants to restart adaptive -> restart.

        Adaptive restart implies picking up the simulation at the last iteration, whereas restart runs PELE from scratch
        but without overwriting adaptive.conf and pele.conf files.
        """
        output_dir = os.path.join(context.parameters.pele_dir, context.parameters.output)

        if context.parameters.adaptive_restart and not os.path.exists(output_dir):
            context.parameters.adaptive_restart = False

    def set_params(self, simulation_type):
        """
        Make sure all simulations are set to False, then set the one you need to True. This is to avoid scenarios where
        we have induced_fit_fast: true and rescoring: true, because some random parameters were carried over.
        """
        for sim in ft.all_simulations:  # make sure everything else is False
            setattr(context.parameters, sim, "")
        setattr(context.parameters, simulation_type, True)  # set the simulation you need

    def set_user_params(self):
        """
        Overriding default pele_env variables by user-defined parameters from input.yaml.
        """
        if self.options:
            for key, value in self.options.items():
                setattr(context.parameters, key, value)

    def create_folders(self):
        context.parameters.create_files_and_folders()

    def water_handler(self):
        """
        In the PPI and Allosteric packages, water perturbation is only executed in the refinement simulation. We
        temporarily hide n_waters parameter to avoid adding water molecules to the global/interface exploration.
        """
        if getattr(context.yaml_parser, "n_waters", None):
            if (self.keyword == "full" and context.parameters_builder.package == "allosteric") or (
                    self.keyword == "induced_fit_exhaustive"
                    and context.parameters_builder.package == "ppi"
            ):
                self.hide_water()
            elif (
                    self.keyword == "induced_fit_exhaustive"
                    and context.parameters_builder.package == "site_finder"
            ) or (self.keyword == "rescoring" and context.parameters_builder.package == "ppi"):
                self.add_water()

    def hide_water(self):
        """
        Temporarily hide n_waters as _n_waters, so that water perturbation is not performed during global/interface
        exploration.
        """
        context.parameters._n_waters = context.parameters.n_waters
        context.parameters.n_waters = 0
        context.parameters.water_arg = None
        context.parameters.waters = None

    def add_water(self):
        """
        Reinstate hidden n_waters.
        """
        context.parameters.n_waters = getattr(context.parameters, "_n_waters", None)
        context.parameters.waters = "all_waters"


class GlobalExploration(Simulation):

    keyword = "full"

    def set_package_params(self):
        pass


class LocalExplorationFast(Simulation):

    keyword = "induced_fit_fast"

    def set_package_params(self):
        pass


class LocalExplorationExhaustive(Simulation):

    keyword = "induced_fit_exhaustive"

    def set_package_params(self):
        if context.parameters_builder.package == "ppi":
            context.parameters.water_arg = None
            context.parameters.system = prepare_structure(
                context.parameters.system,
                context.parameters.ligand_pdb,
                context.parameters.protein,
                remove_water=False,
            )


class Rescoring(Simulation):

    keyword = "rescoring"

    def set_package_params(self):
        if context.parameters_builder.package == "ppi":
            self.set_ppi_params()

    def set_ppi_params(self):
        """
        Overrides default Rescoring parameters if running PPI.
        """
        if not context.parameters.test:
            context.parameters.iterations = 1
            context.parameters.steps = 100
            context.parameters.box_radius = 100

        if hasattr(context.parameters, "_n_waters"):
            context.parameters.n_waters = context.parameters._n_waters


class GPCR(Simulation):

    keyword = "gpcr_orth"

    def set_package_params(self):
        context.parameters.orthosteric_site = context.yaml_parser.orthosteric_site
        context.parameters.initial_site = context.yaml_parser.initial_site
        context.parameters.center_of_interface = context.parameters.initial_site
        box_center, box_radius = helpers.retrieve_box(
            context.parameters.system,
            context.parameters.initial_site,
            context.parameters.orthosteric_site,
            weights=[0.35, 0.65],
        )
        context.parameters.box_center = (
            context.yaml_parser.box_center
            if context.yaml_parser.box_center
            else box_center
        )
        context.parameters.box_radius = (
            context.yaml_parser.box_radius
            if context.yaml_parser.box_radius
            else box_radius
        )
        context.parameters.randomize = True


class OutIn(Simulation):

    keyword = "out_in"

    def set_package_params(self):
        self._check_mandatory_fields()
        self.set_outin_params()

    def _check_mandatory_fields(self):
        compulsory_flags = ["final_site", "initial_site"]
        for flag in compulsory_flags:
            if getattr(context.yaml_parser, flag) is None:
                raise custom_errors.OutInError(
                    f"Flag {flag} must be specified for out_in package."
                )

    def set_outin_params(self):
        context.parameters.final_site = context.yaml_parser.final_site
        context.parameters.initial_site = context.yaml_parser.initial_site
        context.parameters.center_of_interface = context.parameters.initial_site
        box_center, box_radius = helpers.retrieve_box(
            context.yaml_parser.system,
            context.parameters.initial_site,
            context.parameters.final_site,
            weights=[0.35, 0.65],
        )
        context.parameters.box_center = (
            context.yaml_parser.box_center
            if context.yaml_parser.box_center
            else box_center
        )
        context.parameters.box_radius = (
            context.yaml_parser.box_radius
            if context.yaml_parser.box_radius
            else box_radius
        )
        context.parameters.randomize = True


class CovalentDockingExploration(Simulation):

    keyword = "covalent_docking"

    def set_package_params(self):
        context.parameters.residue_type = self.get_residue_type()
        self.correct_system()
        self.set_general_perturbation_params()
        context.parameters.system = self.run_ppp()
        parametrizer.parametrize_covalent_residue(context.parameters.pele_data, context.parameters.pele_dir, context.parameters.gridres,
                                                  context.parameters.residue_type, context.parameters.residue,
                                                  ppp_system=context.parameters.system)

    def set_general_perturbation_params(self):
        """
        Sets parameters for the initial side chain perturbation, making sure we set the correct working folder and
        ignore refinement distance for now.
        """
        if isinstance(context.parameters.skip_ligand_prep, list):
            context.parameters.skip_ligand_prep.append(context.parameters.residue)
        else:
            context.parameters.skip_ligand_prep = [context.parameters.residue]

    def get_residue_type(self):
        """
        Extracts name of the residue the covalent ligand is bound to before correcting the system.
        """
        chain, residue_number = context.parameters.covalent_residue.split(":")
        residue_type = helpers.get_residue_name(context.parameters.system, chain, residue_number)
        return residue_type

    def correct_system(self):
        """
        Moves the covalent ligand to the other residues, replaces HETATM with ATOM, then assigns the resulting PDB
        as the new system. Then adds original CONECT lines back to the extracted LIG.pdb.
        """
        from frag_pele.Covalent import pdb_corrector

        corrected_system = os.path.join(
            context.parameters.inputs_dir,
            os.path.basename(context.parameters.system.replace(".pdb", "_corrected.pdb")),
        )
        chain, residue_number = context.parameters.covalent_residue.split(":")

        pdb_corrector.run(
            context.parameters.system,
            chain,
            int(residue_number),
            corrected_system,
            ligand_resname=context.parameters.residue,
            ligand_chain=context.parameters.chain,
        )

        self.retrieve_ligand_conects()
        context.parameters.system = corrected_system

    def retrieve_ligand_conects(self):
        """
        Maps atom numbers from the original system PDB and adds modified CONECT lines to the extracted ligand PDB.
        """
        original_atom_numbers = list()
        original_conects = list()

        covalent_chain, covalent_resnum = context.parameters.covalent_residue.split(":")
        extracted_ligand = os.path.join(context.parameters.inputs_dir, f"{context.parameters.residue}.pdb")

        # Load original PDB
        with open(context.parameters.system, "r") as system_file:
            system_lines = [
                line
                for line in system_file.readlines()
                if line.startswith("ATOM") or line.startswith("HETATM") or line.startswith("CONECT")]

        # Find ligand covalent residue lines in the original PDB and extract atom numbers
        for line in system_lines:
            if line[17:20].strip() == context.parameters.residue:
                original_atom_numbers.append(line[7:11].strip())
            if line[22:26].strip() == covalent_resnum and line[21] == covalent_chain:
                original_atom_numbers.append(line[7:11].strip())

        # Extract CONECT lines containing relevant atom numbers (residue and ligand)
        for line in system_lines:
            if line.startswith("CONECT") and any(
                    number in line for number in original_atom_numbers
            ):
                original_conects.append(line)

        original_conects = "".join(original_conects)

        # Extract new ligand atom numbers
        with open(extracted_ligand, "r") as ligand_file:
            ligand_lines = ligand_file.readlines()
            new_atom_numbers = [line[7:11].strip() for line in ligand_lines]
            # Schrodinger needs 4 spaces, otherwise it makes a mess
            new_atom_numbers = [f"    {number}" for number in new_atom_numbers]

        # Go over CONECT lines and update atom numbers
        for old, new in zip(original_atom_numbers, new_atom_numbers):
            original_conects = re.sub(rf"\b{old}\b", new, original_conects)

        # Append mapped CONECT lines to the extracted LIG.pdb
        with open(extracted_ligand, "a") as new_ligand:
            new_ligand.write(original_conects)

    def run_ppp(self):
        """
        Retrieves non-standard residue (in this case the covalent ligand) and runs PPP.
        Returns
        -------
            PDB file with preprocessed system.
        """
        import PPP.main as ppp

        context.parameters.nonstandard.extend(helpers.find_nonstd_residue(context.parameters.system))

        if context.parameters.no_ppp:
            return context.parameters.system
        else:
            return ppp.main(
                context.parameters.system,
                context.parameters.inputs_dir,
                output_pdb=["", ],
                charge_terminals=context.parameters.charge_ter,
                no_gaps_ter=context.parameters.gaps_ter,
                mid_chain_nonstd_residue=context.parameters.nonstandard,
                skip=context.parameters.skip_prep,
                back_constr=context.parameters.ca_constr,
                constrain_smiles=None,
                ligand_pdb=context.parameters.ligand_ref,
                ca_interval=context.parameters.ca_interval)[0]


class CovalentDockingRefinement(Simulation):

    keyword = "covalent_docking_refinement"

    def set_package_params(self):
        """
        Sets parameters for the refinement side chain perturbation, including the refinement distance (default = 10 A).
        """
        context.parameters.no_ppp = True
        context.parameters.covalent_docking_refinement = True
        self.recover_templates_from_global()

    def recover_templates_from_global(self):
        """
        Sets templates created in the first part of the simulation as external templates and ensures
        the parametrization of those ligands is skipped during refinement.
        """
        global_pele_dir = os.path.dirname(context.parameters.pele_dir)

        templates = glob.glob(
            os.path.join(global_pele_dir, "*", "DataLocal", "Templates", "OPLS2005", "Protein", "*")
        ) + glob.glob(
            os.path.join(
                global_pele_dir, "DataLocal", "Templates", "OPLS2005", "HeteroAtoms", "*"
            )
        )

        context.parameters.template = [
            template
            for template in templates
            if os.path.basename(template) != "templates_generated"
        ]

        context.parameters.rotamers = glob.glob(
            os.path.join(global_pele_dir, "DataLocal", "LigandRotamerLibs", "*")
        )

        for template in context.parameters.rotamers:
            context.parameters.skip_ligand_prep.append(
                os.path.basename(template.strip(".rot.assign"))
            )

        context.parameters.skip_ligand_prep = list(set(context.parameters.skip_ligand_prep))
