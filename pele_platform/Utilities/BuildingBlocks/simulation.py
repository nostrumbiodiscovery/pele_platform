from abc import abstractmethod
from dataclasses import dataclass
import glob
import os

import pele_platform.Adaptive.simulation as si
import pele_platform.Errors.custom_errors as ce
from pele_platform.Utilities.Helpers.helpers import retrieve_box
from pele_platform.Utilities.BuildingBlocks.preparation import prepare_structure
from pele_platform.Utilities.BuildingBlocks import blocks
import pele_platform.Utilities.Parameters.pele_env as pv
import pele_platform.features.adaptive as ft


@dataclass
class Simulation(blocks.Block):
    """
    Base Simulation class to run all simulation types.
    Both input and output should always be EnviroBuilder.

    One class to rule them all, one class to find them, one class to bring them all and in PELE bind them.
    """

    env: pv.EnviroBuilder
    options: dict
    folder_name: str

    def run(self):
        self.simulation_setup()
        self.env = si.run_adaptive(self.env)
        return self.env

    def simulation_setup(self):
        self.set_params(simulation_type=self.keyword)
        self.set_working_folder()
        self.env.build_adaptive_variables(self.env.initial_args)
        self.set_package_params()
        self.set_user_params()
        self.create_folders()
        if hasattr(self.env, "next_step"):
            self.env.input = glob.glob(self.env.next_step)
        self.restart_checker()
        self.water_handler()

    @abstractmethod
    def set_package_params(self):
        pass

    def restart_checker(self):
        """
        Check if we should run restart based on the presence of the output folder and 'adaptive_restart' flag in
        input.yaml.
        """
        output_dir = os.path.join(self.env.pele_dir, self.env.output)
        if self.env.adaptive_restart and os.path.exists(output_dir):
            self.env.adaptive_restart = True
        else:
            self.env.adaptive_restart = False

    def set_params(self, simulation_type):
        """
        Make sure all simulations are set to False, then set the one you need to True. This is to avoid scenarios where
        we have induced_fit_fast: true and rescoring: true, because some random parameters were carried over.
        """
        for sim in ft.all_simulations:  # make sure everything else is False
            setattr(self.env, sim, False)
        setattr(self.env, simulation_type, True)  # set the simulation you need

    def set_user_params(self):
        """
        Overriding default pele_env variables by user-defined parameters from input.yaml.
        """
        if self.options:
            for key, value in self.options.items():
                setattr(self.env, key, value)

    def set_working_folder(self):
        """
        Check if user specified a custom folder name for this simulation block. If not, use the automatically
        generated one.
        """
        self.original_dir = os.path.abspath(os.getcwd())
        if self.options:
            user_folder = self.options.get("working_folder", None)
            self.env.folder_name = user_folder if user_folder else self.folder_name
        else:
            self.env.folder_name = self.folder_name

        if hasattr(self.env, "pele_dir"):
            self.env.initial_args.folder = os.path.join(
                os.path.dirname(self.env.pele_dir), self.env.folder_name
            )

    def create_folders(self):
        self.env.create_files_and_folders()

    def water_handler(self):
        """
        In the PPI and Allosteric packages, water perturbation is only executed in the refinement simulation. We
        temporarily hide n_waters parameter to avoid adding water molecules to the global/interface exploration.
        """
        if getattr(self.env.initial_args, "n_waters", None):
            if (self.keyword == "full" and self.env.package == "allosteric") or (
                self.keyword == "inducef_fit_exhaustive" and self.env.package == "ppi"
            ):
                self.hide_water()
            elif (
                self.keyword == "induced_fit_exhaustive"
                and self.env.package == "allosteric"
            ) or (self.keyword == "rescoring" and self.env.package == "ppi"):
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


@dataclass
class GlobalExploration(Simulation):
    keyword: str = "full"

    def set_package_params(self):
        pass


@dataclass
class LocalExplorationFast(Simulation):
    keyword: str = "induced_fit_fast"

    def set_package_params(self):
        pass


@dataclass
class LocalExplorationExhaustive(Simulation):
    keyword: str = "induced_fit_exhaustive"

    def set_package_params(self):
        if self.env.package == "ppi":
            self.env.water_arg = None
            self.env.system = prepare_structure(
                self.env.system,
                self.env.ligand_pdb,
                self.env.protein,
                remove_water=False,
            )


@dataclass
class Rescoring(Simulation):
    keyword: str = "rescoring"

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


@dataclass
class GPCR(Simulation):
    keyword: str = "gpcr_orth"

    def set_package_params(self):
        self.env.orthosteric_site = self.env.initial_args.orthosteric_site
        self.env.initial_site = self.env.initial_args.initial_site
        self.env.center_of_interface = self.env.initial_site
        box_center, box_radius = retrieve_box(
            self.env.system,
            self.env.initial_site,
            self.env.orthosteric_site,
            weights=[0.35, 0.65],
        )
        self.env.box_center = (
            self.env.initial_args.box_center
            if self.env.initial_args.box_center
            else box_center
        )
        self.env.box_radius = (
            self.env.initial_args.box_radius
            if self.env.initial_args.box_radius
            else box_radius
        )
        self.env.randomize = True


@dataclass
class OutIn(Simulation):
    keyword: str = "out_in"

    def set_package_params(self):
        self._check_mandatory_fields()
        self.set_outin_params()

    def _check_mandatory_fields(self):
        compulsory_flags = ["final_site", "initial_site"]
        for flag in compulsory_flags:
            if getattr(self.env.initial_args, flag) is None:
                raise ce.OutInError(
                    f"Flag {flag} must be specified for out_in package."
                )

    def set_outin_params(self):
        self.env.final_site = self.env.initial_args.final_site
        self.env.initial_site = self.env.initial_args.initial_site
        self.env.center_of_interface = self.env.initial_site
        box_center, box_radius = retrieve_box(
            self.env.initial_args.system,
            self.env.initial_site,
            self.env.final_site,
            weights=[0.35, 0.65],
        )
        self.env.box_center = (
            self.env.initial_args.box_center
            if self.env.initial_args.box_center
            else box_center
        )
        self.env.box_radius = (
            self.env.initial_args.box_radius
            if self.env.initial_args.box_radius
            else box_radius
        )
        self.env.randomize = True
