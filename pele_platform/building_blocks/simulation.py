from abc import abstractmethod
import glob
import os

import pele_platform.Adaptive.simulation as si
import pele_platform.Errors.custom_errors as ce
from pele_platform.Utilities.Helpers.helpers import retrieve_box
from pele_platform.building_blocks.preparation import prepare_structure
from pele_platform.building_blocks import blocks
import pele_platform.Utilities.Parameters.parameters as pv
import pele_platform.features.adaptive as ft


class Simulation(blocks.Block):
    """
    Base Simulation class to run all simulation types.

    One class to rule them all, one class to find them, one class to bring them all and in PELE bind them.
    """
    def __init__(
        self,
        parameters_builder: pv.ParametersBuilder,
        options: dict,
        folder_name: str,
        env: pv.Parameters = None,
    ) -> None:
        """
        Initialize Simulation class.

        Parameters
        ----------
        parameters_builder : ParametersBuilder
        options : dict
        folder_name : str
        env : Parameters
        """
        self.options = options
        self.folder_name = folder_name
        self.builder = parameters_builder
        self.original_dir = os.path.abspath(os.getcwd())
        self.env = env if env else self.builder.build_adaptive_variables(self.builder.initial_args)

    def run(self) -> (pv.ParametersBuilder, pv.Parameters):
        self.simulation_setup()
        self.env = si.run_adaptive(self.env)
        return self.builder, self.env

    def simulation_setup(self):
        self.set_working_folder()
        self.set_params(simulation_type=self.keyword)
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
            setattr(self.env, sim, "")
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
        if self.options:
            user_folder = self.options.get("working_folder", None)
            self.env.folder_name = user_folder if user_folder else self.folder_name
        else:
            self.env.folder_name = self.folder_name

        if self.env.pele_dir == self.builder.pele_dir:
            self.env.pele_dir = os.path.join(self.builder.pele_dir, self.env.folder_name)
        else:
            self.env.pele_dir = os.path.join(os.path.dirname(self.env.pele_dir), self.env.folder_name)

        self.env.inputs_dir = os.path.join(self.env.pele_dir, "input")

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
                and self.builder.package == "allosteric"
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
    def __init__(
        self, parameters_builder: pv.ParametersBuilder, options: dict, folder_name: str, env: pv.Parameters
    ):
        super().__init__(parameters_builder, options, folder_name, env)
        self.keyword = "full"

    def set_package_params(self):
        pass


class LocalExplorationFast(Simulation):
    def __init__(
        self, parameters_builder: pv.ParametersBuilder, options: dict, folder_name: str, env: pv.Parameters
    ):
        super().__init__(parameters_builder, options, folder_name, env)
        self.keyword = "induced_fit_fast"

    def set_package_params(self):
        pass


class LocalExplorationExhaustive(Simulation):
    def __init__(
        self, parameters_builder: pv.ParametersBuilder, options: dict, folder_name: str, env: pv.Parameters
    ):
        super().__init__(parameters_builder, options, folder_name, env)
        self.keyword = "induced_fit_exhaustive"

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
    def __init__(
        self, parameters_builder: pv.ParametersBuilder, options: dict, folder_name: str, env: pv.Parameters
    ):
        super().__init__(parameters_builder, options, folder_name, env)
        self.keyword = "rescoring"

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
    def __init__(
        self, parameters_builder: pv.ParametersBuilder, options: dict, folder_name: str, env: pv.Parameters
    ):
        super().__init__(parameters_builder, options, folder_name, env)
        self.keyword = "gpcr_orth"

    def set_package_params(self):
        self.env.orthosteric_site = self.builder.initial_args.orthosteric_site
        self.env.initial_site = self.builder.initial_args.initial_site
        self.env.center_of_interface = self.env.initial_site
        box_center, box_radius = retrieve_box(
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
    def __init__(
        self, parameters_builder: pv.ParametersBuilder, options: dict, folder_name: str, env: pv.Parameters
    ):
        super().__init__(parameters_builder, options, folder_name, env)
        self.keyword = "out_in"

    def set_package_params(self):
        self._check_mandatory_fields()
        self.set_outin_params()

    def _check_mandatory_fields(self):
        compulsory_flags = ["final_site", "initial_site"]
        for flag in compulsory_flags:
            if getattr(self.builder.initial_args, flag) is None:
                raise ce.OutInError(
                    f"Flag {flag} must be specified for out_in package."
                )

    def set_outin_params(self):
        self.env.final_site = self.builder.initial_args.final_site
        self.env.initial_site = self.builder.initial_args.initial_site
        self.env.center_of_interface = self.env.initial_site
        box_center, box_radius = retrieve_box(
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
