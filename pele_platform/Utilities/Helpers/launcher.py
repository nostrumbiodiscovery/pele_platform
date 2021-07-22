from dataclasses import dataclass

import pele_platform.DrugDesign.package_launchers
import pele_platform.Checker.main as ck
import pele_platform.Frag.simulation as fr
import pele_platform.Utilities.Parameters.pele_env as pv
from pele_platform.constants import constants
import argparse
import pkg_resources

PACKAGES = dict(
    frag_core=fr.FragRunner,
    ppi=pele_platform.DrugDesign.package_launchers.PPILauncher,
    allosteric=pele_platform.DrugDesign.package_launchers.AllostericLauncher,
    gpcr_orth=pele_platform.DrugDesign.package_launchers.GPCRLauncher,
    out_in=pele_platform.DrugDesign.package_launchers.OutInLauncher,
    adaptive=pele_platform.DrugDesign.package_launchers.AdaptiveLauncher,
    induced_fit_exhaustive=pele_platform.DrugDesign.package_launchers.InducedFitExhaustiveLauncher,
    induced_fit_fast=pele_platform.DrugDesign.package_launchers.InducedFitFastLauncher,
    workflow=pele_platform.DrugDesign.package_launchers.WorkflowLauncher,
)


@dataclass
class Launcher:

    _args: argparse.ArgumentParser
    frag: str = "frag"
    ppi: str = "PPI"
    site_finder: str = "site_finder"
    gpcr_orth: str = "gpcr_orth"
    out_in: str = "out_in"
    adaptive: str = "adaptive"
    saturated_mutagenesis: str = "saturated_mutagenesis"
    interaction_restrictions: str = "interaction_restrictions"
    covalent_docking: str = "covalent_docking"

    def launch(self) -> pv.ParametersBuilder:
        # Launch package from input.yaml
        self.env = pv.ParametersBuilder()
        self.env.initial_args = self._args  # to keep the original args
        print(constants.version_header.format(pkg_resources.get_distribution("pele_platform").version))
        return self.launch_package()

    def launch_package(self) -> pv.EnviroBuilder:
        if not self._args.no_check:
            ck.check_executable_and_env_variables(self._args)
        # Launch package from API
        for package_name, package in PACKAGES.items():
            selected_package = getattr(self._args, package_name)
            if selected_package:
                break
        else:
            package_name = "adaptive"
            package = pele_platform.DrugDesign.package_launchers.AdaptiveLauncher
        self.env.package = self._args.package = package_name
        return package(self.env).run()
