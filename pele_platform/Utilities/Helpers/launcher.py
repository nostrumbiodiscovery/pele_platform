from dataclasses import dataclass
import pele_platform.Checker.main as ck
import pele_platform.Frag.simulation as fr
import pele_platform.Adaptive.simulation as ad
from pele_platform.Allosteric.main import run_allosteric
import pele_platform.gpcr.main as gpcr
import pele_platform.out_in.main as outin
from pele_platform.PPI.main import run_ppi
import pele_platform.Utilities.Parameters.pele_env as pv
import argparse


@dataclass
class Launcher:

    _args: argparse.ArgumentParser
    frag: str="frag"
    ppi: str="PPI"
    allosteric: str="allosteric"
    gpcr_orth: str="gpcr_orth"
    out_in: str="out_in"
    adaptive: str="adaptive"

    def launch(self) -> pv.EnviroBuilder:
        # Launch package from input.yaml
        self._define_package_to_run()
        job_variables = self.launch_package(self._args.package, no_check=self._args.no_check)
        return job_variables

    def launch_package(self, package: str, no_check=False) -> pv.EnviroBuilder:
        # Launch package from API
        if not no_check:
            ck.check_executable_and_env_variables(self._args)
        if package == self.adaptive:
            job_variables = ad.run_adaptive(self._args)
        elif package == self.gpcr_orth:
            job_variables = gpcr.GpcrLauncher(self._args).run_gpcr_simulation()
        elif package == self.out_in:
            job_variables = outin.OutInLauncher(self._args).run_gpcr_simulation()
        elif package == self.allosteric:
            job_variables = run_allosteric(self._args)
        elif package == self.ppi:
            job_variables = run_ppi(self._args)
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
        else: 
            self._args.package = self.adaptive
