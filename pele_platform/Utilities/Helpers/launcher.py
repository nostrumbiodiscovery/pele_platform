from dataclasses import dataclass
import pkg_resources

import pele_platform.drug_design.package_launchers
import pele_platform.enzyme_engineering.saturated_mutagenesis
import pele_platform.Checker.main as checker
from pele_platform.Frag.simulation import FragRunner
from pele_platform.constants import constants
from pele_platform.context import context
from pele_platform.Utilities.Parameters.parameters import ParametersBuilder


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
    covalent_residue=pele_platform.drug_design.package_launchers.CovalentDocking,
)


@dataclass
class Launcher:

    def launch(self):
        """
        Launches PELE package from input.yaml.

        Returns
        -------
            Parameters object(s) with simulation parameters.
        """
        context.parameters_builder = ParametersBuilder()  # move to Context init?
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
        if not context.yaml_parser.no_check:
            checker.check_executable_and_env_variables(context.yaml_parser)

        for package_name, package in PACKAGES.items():
            selected_package = getattr(context.yaml_parser, package_name)
            if selected_package:
                break
        else:
            package_name = "adaptive"
            package = pele_platform.drug_design.package_launchers.AdaptiveLauncher

        context.parameters_builder.package = context.yaml_parser.package = package_name
        return package().run()
