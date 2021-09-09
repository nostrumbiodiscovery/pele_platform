import pkg_resources

import pele_platform.drug_design.pipelines
import pele_platform.enzyme_engineering.saturated_mutagenesis
import pele_platform.Checker.main as checker
from pele_platform.Frag.simulation import FragRunner
from pele_platform.constants import constants
from pele_platform.context import context
from pele_platform.Utilities.Parameters.parameters import ParametersBuilder


PACKAGES = dict(
    frag_core=FragRunner,
    ppi=pele_platform.drug_design.pipelines.PPI,
    site_finder=pele_platform.drug_design.pipelines.SiteFinder,
    gpcr_orth=pele_platform.drug_design.pipelines.GPCR,
    out_in=pele_platform.drug_design.pipelines.OutIn,
    adaptive=pele_platform.drug_design.pipelines.Adaptive,
    induced_fit_exhaustive=pele_platform.drug_design.pipelines.InducedFitExhaustive,
    induced_fit_fast=pele_platform.drug_design.pipelines.InducedFitFast,
    workflow=pele_platform.drug_design.pipelines.Workflow,
    covalent_residue=pele_platform.drug_design.pipelines.CovalentDocking,
)


class Launcher:

    print(
        constants.version_header.format(
            pkg_resources.get_distribution("pele_platform").version
        )
    )

    def launch_pele(self):
        """
        Launches PELE packages from input.yaml

        Returns
        -------
            Parameters object(s) with simulation parameters.
        """
        if not context.yaml_parser.no_check:
            checker.check_executable_and_env_variables(context.yaml_parser)

        for package_name, package in PACKAGES.items():
            selected_package = getattr(context.yaml_parser, package_name)
            if selected_package:
                break
        else:
            package_name = "adaptive"
            package = pele_platform.drug_design.pipelines.Adaptive

        context.parameters_builder = ParametersBuilder()  # move to Context init?
        context.parameters_builder.package = context.yaml_parser.package = package_name

        return package().run()

    def launch_pydock(self):
        pass
