import pele_platform.Checker.main as ck
import pele_platform.Frag.simulation as fr
import pele_platform.Adaptive.simulation as ad
from pele_platform.Allosteric.main import run_allosteric
import pele_platform.gpcr.main as gpcr
from pele_platform.PPI.main import run_ppi

class Launcher():


    def __init__(self, arguments):
        self._args = arguments
        self.frag = "frag"
        self.ppi = "PPI"
        self.allosteric = "allosteric"
        self.gpcr_orth = "gpcr_orth"
        self.adaptive = "adaptive"

    def launch(self):
        self._define_package_to_run()
        if not self._args.no_check:
            ck.check_executable_and_env_variables(self._args)
        if self._args.package == self.adaptive:
            job_variables = ad.run_adaptive(self._args)
        elif self._args.package == self.gpcr_orth:
            job_variables = gpcr.GpcrLauncher(self._args).run_gpcr_simulation()
        elif self._args.package == self.allosteric:
            job_variables = run_allosteric(self._args)
        elif self._args.package == self.ppi:
            job_variables = run_ppi(self._args)
        elif self._args.package == self.frag:
            #Set variables and input ready 
            job_variables = fr.FragRunner(self._args).run_simulation()
        return job_variables


    def _define_package_to_run(self):
        if self._args.frag_core:
            self._args.package = self.frag
        elif self._args.ppi:
            self._args.package = self.ppi
        elif self._args.allosteric:
            self._args.package = self.allosteric
        elif self._args.gpcr_orth:
            self._args.package = self.gpcr_orth
        else: 
            self._args.package = self.adaptive

