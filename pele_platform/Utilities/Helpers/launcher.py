from dataclasses import dataclass
import pele_platform.Utilities.Parameters.pele_env as pele
import pele_platform.Checker.main as ck
import pele_platform.Frag.simulation as fr
import pele_platform.Adaptive.main as adp
import pele_platform.Allosteric.main as al
import pele_platform.GPCR.main as gpcr
import pele_platform.InducedFit.main as ind
import pele_platform.out_in.main as outin
import pele_platform.PPI.main as ppi
import pele_platform.Utilities.Parameters.pele_env as pv
import argparse

PACKAGES = dict(
    frag_core=fr.FragRunner,
    ppi=ppi.PPILauncher,
    allosteric=al.AllostericLauncher,
    gpcr_orth=gpcr.GPCRLauncher,
    out_in=outin.OutInLauncher,
    adaptive=adp.AdaptiveLauncher,
    induced_fit_exhaustive=ind.InducedFitExhaustiveLauncher,
    induced_fit_fast=ind.InducedFitFastLauncher,
    workflow=adp.WorkflowLauncher,
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
            package = adp.AdaptiveLauncher
        self.env.package = self._args.package = package_name
        return package(self.env).run()
