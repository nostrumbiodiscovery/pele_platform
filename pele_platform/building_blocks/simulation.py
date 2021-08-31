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


class Simulation(blocks.Block):
    """
    Base Simulation class to run all simulation types.

    One class to rule them all, one class to find them, one class to bring them all and in PELE bind them.
    """
    keyword = None

    def run(self) -> (parameters.ParametersBuilder, parameters.Parameters):
        self.simulation_setup()
        self.env = si.run_adaptive(self.env)
        self.post_simulation_cleanup()
        return self.builder, self.env

    def simulation_setup(self):
        self.set_working_folder()
        self.set_user_params()
        self.restart_checker()
        self.create_folders()
        self.set_params(simulation_type=self.keyword)
        self.set_package_params()
        if hasattr(self.env, "next_step"):
            self.env.input = glob.glob(self.env.next_step)
        self.water_handler()

    def post_simulation_cleanup(self):
        """
        Sets both types of restart to false, so that the platform does not look for logger file while executing the
        next building block. Disables PDB preprocessing, since it would have been done at the first step.
        """
        self.env.restart = False
        self.env.adaptive_restart = False
        self.env.no_ppp = True

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
        output_dir = os.path.join(self.env.pele_dir, self.env.output)

        if self.env.adaptive_restart and not os.path.exists(output_dir):
            self.env.adaptive_restart = False

    def set_params(self, simulation_type):
        """
        Make sure all simulations are set to False, then set the one you need to True. This is to avoid scenarios where
        we have induced_fit_fast: true and rescoring: true, because some random parameters were carried over.
        """
        for sim in ft.all_simulations:  # make sure everything else is False
            setattr(self.env, sim, "")
        setattr(self.env, simulation_type, True)  # set the simulation you need

    def set_user_params(self):
        """
        Overriding default pele_env variables by user-defined parameters from input.yaml.
        """
        if self.options:
            for key, value in self.options.items():
                setattr(self.env, key, value)

    def create_folders(self):
        self.env.create_files_and_folders()

    def water_handler(self):
        """
        In the PPI and Allosteric packages, water perturbation is only executed in the refinement simulation. We
        temporarily hide n_waters parameter to avoid adding water molecules to the global/interface exploration.
        """
        if getattr(self.builder.initial_args, "n_waters", None):
            if (self.keyword == "full" and self.builder.package == "allosteric") or (
                    self.keyword == "induced_fit_exhaustive"
                    and self.builder.package == "ppi"
            ):
                self.hide_water()
            elif (
                    self.keyword == "induced_fit_exhaustive"
                    and self.builder.package == "site_finder"
            ) or (self.keyword == "rescoring" and self.builder.package == "ppi"):
                self.add_water()

    def hide_water(self):
        """
        Temporarily hide n_waters as _n_waters, so that water perturbation is not performed during global/interface
        exploration.
        """
        self.env._n_waters = self.env.n_waters
        self.env.n_waters = 0
        self.env.water_arg = None
        self.env.waters = None

    def add_water(self):
        """
        Reinstate hidden n_waters.
        """
        self.env.n_waters = getattr(self.env, "_n_waters", None)
        self.env.waters = "all_waters"


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
        if self.env.package == "ppi":
            self.env.water_arg = None
            self.env.system = prepare_structure(
                self.env.system,
                self.env.ligand_pdb,
                self.env.protein,
                remove_water=False,
            )


class Rescoring(Simulation):

    keyword = "rescoring"

    def set_package_params(self):
        if self.env.package == "ppi":
            self.set_ppi_params()

    def set_ppi_params(self):
        """
        Overrides default Rescoring parameters if running PPI.
        """
        if not self.env.test:
            self.env.iterations = 1
            self.env.steps = 100
            self.env.box_radius = 100

        if hasattr(self.env, "_n_waters"):
            self.env.n_waters = self.env._n_waters


class GPCR(Simulation):

    keyword = "gpcr_orth"

    def set_package_params(self):
        self.env.orthosteric_site = self.builder.initial_args.orthosteric_site
        self.env.initial_site = self.builder.initial_args.initial_site
        self.env.center_of_interface = self.env.initial_site
        box_center, box_radius = helpers.retrieve_box(
            self.env.system,
            self.env.initial_site,
            self.env.orthosteric_site,
            weights=[0.35, 0.65],
        )
        self.env.box_center = (
            self.builder.initial_args.box_center
            if self.builder.initial_args.box_center
            else box_center
        )
        self.env.box_radius = (
            self.builder.initial_args.box_radius
            if self.builder.initial_args.box_radius
            else box_radius
        )
        self.env.randomize = True


class OutIn(Simulation):

    keyword = "out_in"

    def set_package_params(self):
        self._check_mandatory_fields()
        self.set_outin_params()

    def _check_mandatory_fields(self):
        compulsory_flags = ["final_site", "initial_site"]
        for flag in compulsory_flags:
            if getattr(self.builder.initial_args, flag) is None:
                raise custom_errors.OutInError(
                    f"Flag {flag} must be specified for out_in package."
                )

    def set_outin_params(self):
        self.env.final_site = self.builder.initial_args.final_site
        self.env.initial_site = self.builder.initial_args.initial_site
        self.env.center_of_interface = self.env.initial_site
        box_center, box_radius = helpers.retrieve_box(
            self.builder.initial_args.system,
            self.env.initial_site,
            self.env.final_site,
            weights=[0.35, 0.65],
        )
        self.env.box_center = (
            self.builder.initial_args.box_center
            if self.builder.initial_args.box_center
            else box_center
        )
        self.env.box_radius = (
            self.builder.initial_args.box_radius
            if self.builder.initial_args.box_radius
            else box_radius
        )
        self.env.randomize = True


class CovalentDockingExploration(Simulation):

    keyword = "covalent_docking"

    def set_package_params(self):
        self.env.residue_type = self.get_residue_type()
        self.correct_system()
        self.set_general_perturbation_params()
        self.env.system = self.run_ppp()
        parametrizer.parametrize_covalent_residue(self.env.pele_data, self.env.pele_dir, self.env.gridres,
                                                  self.env.residue_type, self.env.residue,
                                                  ppp_system=self.env.system)

    def set_general_perturbation_params(self):
        """
        Sets parameters for the initial side chain perturbation, making sure we set the correct working folder and
        ignore refinement distance for now.
        """
        if isinstance(self.env.skip_ligand_prep, list):
            self.env.skip_ligand_prep.append(self.env.residue)
        else:
            self.env.skip_ligand_prep = [self.env.residue]

    def get_residue_type(self):
        """
        Extracts name of the residue the covalent ligand is bound to before correcting the system.
        """
        chain, residue_number = self.env.covalent_residue.split(":")
        residue_type = helpers.get_residue_name(self.env.system, chain, residue_number)
        return residue_type

    def correct_system(self):
        """
        Moves the covalent ligand to the other residues, replaces HETATM with ATOM, then assigns the resulting PDB
        as the new system. Then adds original CONECT lines back to the extracted LIG.pdb.
        """
        from frag_pele.Covalent import pdb_corrector

        corrected_system = os.path.join(
            self.env.inputs_dir,
            os.path.basename(self.env.system.replace(".pdb", "_corrected.pdb")),
        )
        chain, residue_number = self.env.covalent_residue.split(":")

        pdb_corrector.run(
            self.env.system,
            chain,
            int(residue_number),
            corrected_system,
            ligand_resname=self.env.residue,
            ligand_chain=self.env.chain,
        )

        self.retrieve_ligand_conects()
        self.env.system = corrected_system

    def retrieve_ligand_conects(self):
        """
        Maps atom numbers from the original system PDB and adds modified CONECT lines to the extracted ligand PDB.
        """
        original_atom_numbers = list()
        original_conects = list()

        covalent_chain, covalent_resnum = self.env.covalent_residue.split(":")
        extracted_ligand = os.path.join(self.env.inputs_dir, f"{self.env.residue}.pdb")

        # Load original PDB
        with open(self.env.system, "r") as system_file:
            system_lines = [
                line
                for line in system_file.readlines()
                if line.startswith("ATOM") or line.startswith("HETATM") or line.startswith("CONECT")]

        # Find ligand covalent residue lines in the original PDB and extract atom numbers
        for line in system_lines:
            if line[17:20].strip() == self.env.residue:
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

        self.env.nonstandard.extend(helpers.find_nonstd_residue(self.env.system))

        if self.env.no_ppp:
            return self.env.system
        else:
            return ppp.main(
                self.env.system,
                self.env.inputs_dir,
                output_pdb=["", ],
                charge_terminals=self.env.charge_ter,
                no_gaps_ter=self.env.gaps_ter,
                mid_chain_nonstd_residue=self.env.nonstandard,
                skip=self.env.skip_prep,
                back_constr=self.env.ca_constr,
                constrain_smiles=None,
                ligand_pdb=self.env.ligand_ref,
                ca_interval=self.env.ca_interval)[0]


class CovalentDockingRefinement(Simulation):

    keyword = "covalent_docking_refinement"

    def set_package_params(self):
        """
        Sets parameters for the refinement side chain perturbation, including the refinement distance (default = 10 A).
        """
        self.env.no_ppp = True
        self.env.covalent_docking_refinement = True
        self.recover_templates_from_global()

    def recover_templates_from_global(self):
        """
        Sets templates created in the first part of the simulation as external templates and ensures
        the parametrization of those ligands is skipped during refinement.
        """
        global_pele_dir = os.path.dirname(self.env.pele_dir)

        templates = glob.glob(
            os.path.join(global_pele_dir, "*", "DataLocal", "Templates", "OPLS2005", "Protein", "*")
        ) + glob.glob(
            os.path.join(
                global_pele_dir, "DataLocal", "Templates", "OPLS2005", "HeteroAtoms", "*"
            )
        )

        self.env.template = [
            template
            for template in templates
            if os.path.basename(template) != "templates_generated"
        ]

        self.env.rotamers = glob.glob(
            os.path.join(global_pele_dir, "DataLocal", "LigandRotamerLibs", "*")
        )

        for template in self.env.rotamers:
            self.env.skip_ligand_prep.append(
                os.path.basename(template.strip(".rot.assign"))
            )

        self.env.skip_ligand_prep = list(set(self.env.skip_ligand_prep))
