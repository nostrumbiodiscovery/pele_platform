from dataclasses import dataclass
import pele_platform.Utilities.Parameters.pele_env as pele
import pele_platform.Checker.main as ck
import pele_platform.Frag.simulation as fr
import pele_platform.Adaptive.simulation as ad
import pele_platform.Allosteric.main as al
import pele_platform.GPCR.main as gpcr
import pele_platform.InducedFit.main as ind
import pele_platform.out_in.main as outin
import pele_platform.PPI.main as ppi
import pele_platform.Utilities.Parameters.pele_env as pv
import argparse


@dataclass
class Launcher:

    _args: argparse.ArgumentParser
    frag: str = "frag"
    ppi: str = "PPI"
    allosteric: str = "allosteric"
    gpcr_orth: str = "gpcr_orth"
    out_in: str = "out_in"
    adaptive: str = "adaptive"
    induced_fit_exhaustive: str = "induced_fit_exhaustive"
    induced_fit_fast: str = "induced_fit_fast"


    def launch(self) -> pv.EnviroBuilder:
        # Launch package from input.yaml
        self._define_package_to_run()
        self.env = pele.EnviroBuilder()
        self.env.initial_args = self._args  # to keep the original args
        job_variables = self.launch_package(self._args.package, no_check=self._args.no_check)
        return job_variables

    def launch_package(self, package: str, no_check=False) -> pv.EnviroBuilder:
        # Launch package from API
        if not no_check:
            ck.check_executable_and_env_variables(self._args)
        if package == self.adaptive:
            job_variables = ad.run_adaptive(self._args)
        elif package == self.gpcr_orth:
            job_variables = gpcr.GPCRLauncher(self.env).run()
        elif package == self.out_in:
            job_variables = outin.OutInLauncher(self.env).run()
        elif package == self.allosteric:
            job_variables = al.AllostericLauncher(self.env).run()
        elif package == self.ppi:
            job_variables = ppi.run(self._args)
        elif package == self.induced_fit_fast:
            job_variables = ind.InducedFitFastLauncher(self.env).run()
        elif package == self.induced_fit_exhaustive:
            job_variables = ind.InducedFitExhaustiveLauncher(self.env).run()
        elif package == self.frag:
            # Set variables and input ready
            job_variables = fr.FragRunner(self._args).run_simulation()
        return job_variables

    def _define_package_to_run(self) -> None:
        # Define package being run from input.yaml flags
        if self._args.frag_core:
            self._args.package = self.frag
        elif self._args.ppi:
            self._args.package = self.ppi
        elif self._args.allosteric:
            self._args.package = self.allosteric
        elif self._args.gpcr_orth:
            self._args.package = self.gpcr_orth
        elif self._args.out_in:
            self._args.package = self.out_in
        elif self._args.induced_fit_fast:
            self._args.package = self.induced_fit_fast
        elif self._args.induced_fit_exhaustive:
            self._args.package = self.induced_fit_exhaustive
        else: 
            self._args.package = self.adaptive
