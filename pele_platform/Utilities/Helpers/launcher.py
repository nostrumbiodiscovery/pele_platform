from dataclasses import dataclass

import pele_platform.DrugDesign.package_launchers
import pele_platform.Utilities.Parameters.pele_env as pele
import pele_platform.Checker.main as ck
import pele_platform.Frag.simulation as fr
import pele_platform.Utilities.Parameters.pele_env as pv
import argparse

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

    def launch(self) -> pv.EnviroBuilder:
        # Launch package from input.yaml
        self.env = pele.EnviroBuilder()
        self.env.initial_args = self._args  # to keep the original args
        return self.launch_package()

    def launch_package(self) -> pv.EnviroBuilder:
        if not self._args.no_check:
            ck.check_executable_and_env_variables(self._args)
        # Launch package from API
        for package_name, package in PACKAGES.items():
            if getattr(self._args, package_name):
                break
        else:
            package_name = "adaptive"
            package = pele_platform.DrugDesign.package_launchers.AdaptiveLauncher
        self.env.package = self._args.package = package_name
        return package(self.env).run()
