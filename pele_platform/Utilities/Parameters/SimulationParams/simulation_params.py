import random
import os
import pele_platform.constants.constants as cs
from pele_platform.Utilities.Parameters.SimulationParams.MSMParams import msm_params
from pele_platform.Utilities.Parameters.SimulationParams.GlideParams import glide_params
from pele_platform.Utilities.Parameters.SimulationParams.BiasParams import bias_params
from pele_platform.Utilities.Parameters.SimulationParams.InOutParams import inout_params
from pele_platform.Utilities.Parameters.SimulationParams.WaterExp import waterexp_params
from pele_platform.Utilities.Parameters.SimulationParams.PCA import pca
import pele_platform.Utilities.Helpers.helpers as hp

LOGFILE = '"simulationLogPath" : "$OUTPUT_PATH/logFile.txt",'


class SimulationParams(msm_params.MSMParams, glide_params.GlideParams, bias_params.BiasParams, 
    inout_params.InOutParams,  waterexp_params.WaterExp, pca.PCAParams):


    def __init__(self, args):
        self.simulation_type(args)
        self.main_pele_params(args)
        self.main_adaptive_params(args)
        self.optative_params(args)
        self.system_preparation_params(args)
        self.ligand_params(args)
        self.water_params(args)
        self.box_params(args)
        self.metrics_params(args)
        self.output_params(args)
        self.analysis_params(args)

        #Create all simulation types (could be more efficient --> chnage in future) 
        msm_params.MSMParams.__init__(self, args)
        glide_params.GlideParams.__init__(self, args)
        bias_params.BiasParams.__init__(self, args)
        inout_params.InOutParams.__init__(self, args)
        waterexp_params.WaterExp.__init__(self, args)
        pca.PCAParams.__init__(self, args)


    def simulation_type(self, args):
        self.pele = args.pele if args.adaptive else None
        self.adaptive = args.adaptive if args.adaptive else None
        self.hbond_donor, self.hbond_acceptor = args.hbond

    def main_pele_params(self,args):
        self.system = args.system
        self.residue = args.residue
        self.chain = args.chain
        self.debug = True if args.debug else False
        self.pele_steps = args.pele_steps if args.pele_steps else self.simulation_params.get("pele_steps", None)
        self.logfile =  LOGFILE if args.log else ""
        self.license = '''"{}"'''.format(cs.LICENSE)
        self.anm_freq = self.simulation_params.get("anm_freq", args.anm_freq)
        self.anm_displacement = self.simulation_params.get("anm_displacement", args.anm_displacement)
        self.anm_modes_change = self.simulation_params.get("anm_modes_change", args.anm_modes_change)
        self.sidechain_freq = self.simulation_params.get("sidechain_freq", args.sidechain_freq)
        self.min_freq = self.simulation_params.get("min_freq", args.min_freq)
        self.water_freq = args.water_freq
        self.temperature = self.simulation_params.get("temperature", args.temperature)
        self.sidechain_resolution = args.sidechain_resolution
        self.proximityDetection = "false" if not args.proximityDetection else "true"
        self.steric_trials = args.steric_trials if args.steric_trials else self.simulation_params.get("steric_trials", None)
        self.ca_constr = args.ca_constr
        self.overlap_factor = args.overlap_factor if args.overlap_factor else self.simulation_params.get("overlap_factor", None)
        self.perturbation = args.perturbation if args.perturbation else ""
        self.perturbation_params(args)


    def perturbation_params(self, args):
        if self.perturbation:
            self.selection_to_perturb = args.selection_to_perturb if args.selection_to_perturb else ""
            self.parameters = self.simulation_params.get("params", "") if args.parameters else ""
            self.ligand = cs.LIGAND if self.perturbation else ""
            self.binding_energy = args.binding_energy if args.binding_energy else ""
            self.sasa = args.sasa if args.sasa else ""
        else:
            self.selection_to_perturb = ""
            self.parameters = ""
            self.ligand = ""
            self.binding_energy = ""
            self.sasa = ""


    def main_adaptive_params(self, args):
        self.spawning = args.spawning if args.spawning else self.simulation_params.get("spawning_type", None)
        self.density = args.density if args.density else self.simulation_params.get("density", None)
        self.simulation_type = args.simulation_type if args.simulation_type else self.simulation_params.get("simulation_type", None)
        self.iterations = args.iterations if args.iterations else self.simulation_params.get("iterations", None)
        self.cluster_values = args.cluster_values if args.cluster_values else self.simulation_params.get("cluster_values", None)
        self.cluster_conditions = args.cluster_conditions if args.cluster_conditions else self.simulation_params.get("cluster_conditions", None)
        self.seed = args.seed if args.seed else random.randrange(1, 70000)
        self.templates = os.path.abspath(os.path.join(os.path.dirname(os.path.dirname(__file__)), "PeleTemplates"))
        self.usesrun = "true" if args.usesrun else "false"

    def optative_params(self, args):
        self.input = args.input
        self.forcefield = args.forcefield
        self.solvent = args.solvent
        self.verbose = "true" if args.verbose else "false"
        self.cpus = args.cpus = args.cpus 
        self.restart = args.restart
        self.test = args.test
        self.equil_steps = 1 if self.test else int(args.eq_steps/self.cpus) + 1 #+1 to avoid being 0
        self.equilibration = "true" if args.equilibration else "false"
        self.adaptive_restart = args.adaptive_restart
        self.poses = args.poses
        self.pele_exec = args.pele_exec if args.pele_exec else cs.PELE_BIN
        self.pele_data = args.pele_data if args.pele_data else os.path.join(cs.PELE, "Data")
        self.pele_documents = args.pele_documents if args.pele_documents else os.path.join(cs.PELE, "Documents")

    def system_preparation_params(self, args):
        self.skip_prep = args.skip_prep
        self.nonstandard = args.nonstandard
        self.constraints = None
        self.constrain_smiles = args.constrain_smiles
        self.no_ppp = args.no_ppp

    def ligand_params(self, args):
        self.mae_lig = args.mae_lig
        self.external_template = args.template
        self.external_rotamers = args.rotamers
        self.skip_ligand_prep = args.skip_ligand_prep

    def water_params(self, args):
        self.water_temp = args.water_temp
        self.water_overlap = args.water_overlap
        self.water_constr = args.water_constr
        self.water_trials = args.water_trials
        if args.water_lig or args.water_exp:
            water_arg = args.water_lig if args.water_lig else args.water_exp
            if "all_waters" in [args.water_lig, args.water_exp]:
                water_arg = hp.retrieve_all_waters(self.system)
            self.water_energy = "\n".join([ cs.WATER_ENERGY.format(water.split(":")[0]) for water in water_arg ])
            self.water_energy = None
            self.water = ",".join(['"'+water+'"' for water in water_arg])
            self.water_radius = 6
            # If there is no given center look for it
            if args.water_center:
                self.water_center =  ("[" + ",".join([str(coord) for coord in args.water_center]) + "]")
            else:
                cms = [ hp.find_coords(self.system, water.split(":")[1], water.split(":")[0]) for water in water_arg]
                try:
                    cm = [coord for coord in hp.find_centroid(cms)]
                except TypeError:
                    raise TypeError("Check the specified waters exist")
                self.water_center = cm
                self.water_radius = 6 if  self.water else None
            self.waters = ",".join([ '"' + water + '"' for water in water_arg] )
            self.water = cs.WATER.format(self.water_radius, self.water_center, self.waters, self.water_temp, 
            self.water_trials, self.water_overlap, self.water_constr)
        else:
            self.water_energy = None
            self.water = None
            self.water_radius = None
            self.water_center = None
            self.water = ""


    def box_params(self, args):
        self.box_radius = args.box_radius if args.box_radius else self.simulation_params.get("box_radius", None)
        self.box_center = "["+ ",".join([str(coord) for coord in args.user_center]) + "]" if args.user_center else None

    def metrics_params(self, args):
        self.metrics = None
        self.native = cs.NATIVE.format(os.path.abspath(args.native), self.chain) if args.native else ""
        self.atom_dist = args.atom_dist
        #self.atom_dist = [[args.atom_dist[i], args.atom_dist[i+1]] for i, atom in enumerate(args.atom_dist) if i % 2 == 0] if args.atom_dist else None
        #self.tags = ["distance_{}_to_{}".format(atom1, atom2) for (atom1, atom2) in self.atom_dist] if self.atom_dist else None

    def output_params(self, args):
        self.folder = args.folder
        self.output = args.output
        self.report_name = args.report_name
        self.traj_name = args.traj_name
        self.xtc = args.traj_name.endswith(".xtc")
        self.pdb = args.traj_name.endswith(".pdb")

    def analysis_params(self, args):
        self.analyse = args.analyse
        self.mae = args.mae

