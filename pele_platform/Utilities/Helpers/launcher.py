from dataclasses import dataclass
import pkg_resources

import pele_platform.drug_design.package_launchers
import pele_platform.enzyme_engineering.saturated_mutagenesis
import pele_platform.Checker.main as checker
from pele_platform.Frag.simulation import FragRunner
from pele_platform.Utilities.Parameters.parameters import ParametersBuilder
from pele_platform.constants import constants
from pele_platform.Utilities.Helpers import yaml_parser


PACKAGES = dict(
    frag_core=FragRunner,
    ppi=pele_platform.drug_design.package_launchers.PPILauncher,
    site_finder=pele_platform.drug_design.package_launchers.SiteFinderLauncher,
    gpcr_orth=pele_platform.drug_design.package_launchers.GPCRLauncher,
    out_in=pele_platform.drug_design.package_launchers.OutInLauncher,
    adaptive=pele_platform.drug_design.package_launchers.AdaptiveLauncher,
    induced_fit_exhaustive=pele_platform.drug_design.package_launchers.InducedFitExhaustiveLauncher,
    induced_fit_fast=pele_platform.drug_design.package_launchers.InducedFitFastLauncher,
    workflow=pele_platform.drug_design.package_launchers.WorkflowLauncher,
    saturated_mutagenesis=pele_platform.enzyme_engineering.saturated_mutagenesis.SaturatedMutagenesis,
)


@dataclass
class Launcher:

    _args: yaml_parser.YamlParser
    parameters: ParametersBuilder = None

    def launch(self):
        """
        Launches PELE package from input.yaml.

        Returns
        -------
            ParametersBuilder object(s) with simulation parameters.
        """
        self.parameters = ParametersBuilder()
        self.parameters.initial_args = self._args  # to keep the original args
        print(
            constants.version_header.format(
                pkg_resources.get_distribution("pele_platform").version
            )
        )
        return self.launch_package()

    def launch_package(self):
        """
        Checks which package to run based on keywords in YAML and launches it.

        Returns
        -------
            Parameters object(s) containing simulation parameters for each block.
        """
        if not self._args.no_check:
            checker.check_executable_and_env_variables(self._args)

        for package_name, package in PACKAGES.items():
            selected_package = getattr(self._args, package_name)
            if selected_package:
                break
        else:
            package_name = "adaptive"
            package = pele_platform.drug_design.package_launchers.AdaptiveLauncher

        self.parameters.package = self._args.package = package_name
        return package(self.parameters).run()
