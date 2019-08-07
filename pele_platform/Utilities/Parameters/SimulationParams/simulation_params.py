import random
import os
import pele_platform.constants as cs
from pele_platform.Utilities.Parameters.SimulationParams.MSMParams import msm_params
from pele_platform.Utilities.Parameters.SimulationParams.GlideParams import glide_params
from pele_platform.Utilities.Parameters.SimulationParams.BiasParams import bias_params
import pele_platform.Utilities.Helpers.helpers as hp


class SimulationParams(msm_params.MSMParams, glide_params.GlideParams, bias_params.BiasParams):


    def __init__(self, args):
        self.simulation_type(args)
        self.main_params(args)
        self.optative_params(args)
        self.system_preparation_params(args)
        self.ligand_params(args)
        self.water_params(args)
        self.box_params(args)
        self.metrics_params(args)
        self.output_params(args)

        #Create all simulation types (could be more efficient --> chnage in future) 
        msm_params.MSMParams.__init__(self, args)
        glide_params.GlideParams.__init__(self, args)
        bias_params.BiasParams.__init__(self, args)


    def simulation_type(self, args):
        self.software = args.software
        self.pele = args.pele if args.adaptive else None
        self.adaptive = args.adaptive if args.adaptive else None
        self.hbond_donor, self.hbond_acceptor = args.hbond


    def main_params(self, args):
        self.system = args.system
        self.residue = args.residue
        self.chain = args.chain
        self.adapt_conf = args.adapt_conf
        self.confile = args.confile
        self.seed = random.randrange(1, 70000)
        self.license = '''"{}"'''.format(cs.LICENSE)
        self.templates = os.path.abspath(os.path.join(os.path.dirname(os.path.dirname(__file__)), "PeleTemplates"))

    def optative_params(self, args):
        self.input = args.input
        self.forcefield = args.forcefield
        self.iterations = args.iterations
        self.pele_steps = args.pele_steps
        self.solvent = args.solvent
        self.cpus = args.cpus = args.cpus if not args.test else 4
        self.restart = args.restart
        self.test = args.test
        self.equil_steps = 1 if self.test else int(cs.EQ_STEPS/self.cpus) + 1 #+1 to avoid being 0
        self.equilibration = "true" if args.noeq else "false"

    def system_preparation_params(self, args):
        self.skip_prep = args.skip_prep
        self.nonstandard = args.nonstandard
        self.constraints = None
        self.no_ppp = args.no_ppp

    def ligand_params(self, args):
        self.mae_lig = args.mae_lig
        self.external_template = args.template
        self.external_rotamers = args.rotamers

    def water_params(self, args):
        if args.water_exp:
            self.water_energy = args.water_exp[0].split(":")[0]
            self.water = args.water_exp
        elif args.water_lig:
            for water in args.water_lig:
                print(water.split(":")[0])
            self.water_energy = "\n".join([ cs.WATER_ENERGY.format(water.split(":")[0]) for water in args.water_lig ])
            self.water = ",".join(['"'+water+'"' for water in args.water_lig])
        else:
            self.water_energy = None
            self.water = None
        self.water_radius = 5 if  self.water else None
        self.water_center =  ("[" + ",".join([coord for coord in args.water_center]) + "]") if args.water_center else None

    def box_params(self, args):
        self.box_radius = None
        self.box_center = "["+ ",".join(args.user_center) + "]" if args.user_center else None

    def metrics_params(self, args):
        self.metrics = None
        self.native = cs.NATIVE.format(os.path.abspath(args.native), self.chain) if args.native else ""
        self.atom_dist = args.atom_dist

    def output_params(self, args):
        self.folder = args.folder
        self.report_name = args.report_name
        self.traj_name = args.traj_name
        self.xtc = args.traj_name.endswith(".xtc")
        self.pdb = args.traj_name.endswith(".pdb")

