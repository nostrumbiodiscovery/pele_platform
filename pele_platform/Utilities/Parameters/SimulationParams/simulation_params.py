import random
import os
import glob
import pele_platform.constants.constants as cs
from pele_platform.Utilities.Parameters.SimulationParams.MSMParams import msm_params
from pele_platform.Utilities.Parameters.SimulationParams.GlideParams import glide_params
from pele_platform.Utilities.Parameters.SimulationParams.BiasParams import bias_params
from pele_platform.Utilities.Parameters.SimulationParams.InOutParams import inout_params
from pele_platform.Utilities.Parameters.SimulationParams.WaterExp import waterexp_params
from pele_platform.Utilities.Parameters.SimulationParams.PCA import pca
from pele_platform.Utilities.Parameters.SimulationParams.PPI import ppi
from pele_platform.Utilities.Parameters.SimulationParams.RNA import rna
import pele_platform.Utilities.Helpers.helpers as hp
import pele_platform.Utilities.Helpers.metrics as mt

LOGFILE = '"simulationLogPath" : "$OUTPUT_PATH/logFile.txt",'


class SimulationParams(msm_params.MSMParams, glide_params.GlideParams, bias_params.BiasParams, 
    inout_params.InOutParams,  waterexp_params.WaterExp, pca.PCAParams, ppi.PPIParams, rna.RNAParams):


    def __init__(self, args):
        self.simulation_type(args)
        self.main_pele_params(args)
        self.main_adaptive_params(args)
        self.optative_params(args)
        self.system_preparation_params(args)
        self.ligand_params(args)
        self.anm_params(args)
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
        ppi.PPIParams.__init__(self, args)
        rna.RNAParams.__init__(self, args)


    def simulation_type(self, args):
        self.adaptive = True if args.pele_feature in ["PPI", "adaptive"]  else None
        self.frag_pele = True if args.pele_feature == "frag" else None
        # Trick to let frag handle control fodler parameters --> Improve
        self.complexes = "$PDB" if self.software == "Frag" else "$COMPLEXES"
        self.frag_pele_steps = "$STEPS" if self.software == "Frag" else "$PELE_STEPS"
        self.output_path = "$RESULTS_PATH" if self.software == "Frag" else "$OUTPUT_PATH"

    def main_pele_params(self,args):
        if "*" in args.system:
            self.system = glob.glob(args.system)[0]
            args.input = glob.glob(args.system)
        else:
            self.system = args.system
        self.residue = args.residue
        self.chain = args.chain
        if self.adaptive:
            assert self.system and self.residue and self.chain, "User must define input, residue and chain"
        self.debug = args.debug if args.debug else False
        self.pele_steps = args.pele_steps if args.pele_steps else self.simulation_params.get("pele_steps", 8)
        self.logfile =  LOGFILE if args.log else ""
        self.license = args.pele_license if args.pele_license else os.path.join(cs.PELE, "licenses")
        self.anm_freq = args.anm_freq if args.anm_freq is not None else self.simulation_params.get("anm_freq", 4)
        self.sidechain_freq = args.sidechain_freq if args.sidechain_freq is not None else self.simulation_params.get("sidechain_freq", 2)
        self.min_freq = args.min_freq if args.min_freq is not None else self.simulation_params.get("min_freq", 1)
        self.water_freq = args.water_freq if args.water_freq is not None else self.simulation_params.get("water_freq", 1)
        self.temperature = args.temperature if args.temperature else self.simulation_params.get("temperature", 1500)
        self.sidechain_resolution = args.sidechain_resolution if args.sidechain_resolution else self.simulation_params.get("sidechain_resolution", 30)
        self.proximityDetection = "false" if args.proximityDetection is False else self.simulation_params.get("proximityDetection", "true")
        self.steric_trials = args.steric_trials if args.steric_trials else self.simulation_params.get("steric_trials", 250)
        self.ca_constr = args.ca_constr if args.ca_constr is not None else self.simulation_params.get("ca_constr", 0.5)
        self.overlap_factor = args.overlap_factor if args.overlap_factor else self.simulation_params.get("overlap_factor", 0.65)
        self.perturbation = "" if args.perturbation is False else self.simulation_params.get("perturbation", cs.PERTURBATION)
        self.perturbation_params(args)
        self.com = args.com if args.com else self.simulation_params.get("COMligandConstraint", 0)

    def anm_params(self, args):
        self.anm_displacement = args.anm_displacement if args.anm_displacement else self.simulation_params.get("anm_displacement", 0.75)
        self.anm_modes_change = args.anm_modes_change if args.anm_modes_change else self.simulation_params.get("anm_modes_change", 4)
        self.anm_direction = args.anm_direction if args.anm_direction else self.simulation_params.get("anm_direction", "random")
        self.anm_mix_modes = args.anm_mix_modes if args.anm_mix_modes else self.simulation_params.get("anm_mix_modes", "mixMainModeWithOthersModes")
        self.anm_picking_mode = args.anm_picking_mode if args.anm_picking_mode else self.simulation_params.get("anm_picking_mode", "RANDOM_MODE")
        self.anm_num_of_modes = args.anm_num_of_modes  if args.anm_num_of_modes else self.simulation_params.get("anm_num_of_modes", 6)
        self.anm_relaxation_constr = args.anm_relaxation_constr if args.anm_relaxation_constr else self.simulation_params.get("anm_relaxation_constr", 0.5)
        self.remove_constraints = args.remove_constraints if args.remove_constraints is not None else self.simulation_params.get("remove_constraints", False)

    def perturbation_params(self, args):
        if self.perturbation:
            self.selection_to_perturb = args.selection_to_perturb if args.selection_to_perturb else self.simulation_params.get("selection_to_perturb", cs.SELECTION_TO_PERTURB)
            self.parameters = args.parameters  if args.parameters else self.simulation_params.get("params", True)
            self.ligand = cs.LIGAND if self.perturbation else ""
            self.binding_energy = args.binding_energy if args.binding_energy else self.simulation_params.get("binding_energy", cs.BE)
            self.sasa = args.sasa if args.sasa else self.simulation_params.get("sasa", cs.SASA)
        else:
            self.selection_to_perturb = ""
            self.parameters = ""
            self.ligand = ""
            self.binding_energy = ""
            self.sasa = ""


    def main_adaptive_params(self, args):
        self.spawning = args.spawning if args.spawning else self.simulation_params.get("spawning_type", "independent")
        self.spawning_condition = args.spawning_condition if args.spawning_condition else self.simulation_params.get("spawning_condition", None)
        self.spawning_condition = '"condition": "{}",'.format(self.spawning_condition) if self.spawning_condition else ""
        self.density = args.density if args.density else self.simulation_params.get("density", "null")
        self.simulation_type = args.simulation_type if args.simulation_type else self.simulation_params.get("simulation_type", "pele")
        iterations = 1 if self.spawning == "independent" else 30
        self.iterations = args.iterations if args.iterations else self.simulation_params.get("iterations", iterations)
        self.cluster_values = args.cluster_values if args.cluster_values else self.simulation_params.get("cluster_values", "[1.75, 2.5, 4, 6]")
        self.cluster_conditions = args.cluster_conditions if args.cluster_conditions else self.simulation_params.get("cluster_conditions", "[1, 0.6, 0.4, 0.0]")
        self.seed = args.seed if args.seed else random.randrange(1, 70000)
        self.templates = os.path.abspath(os.path.join(os.path.dirname(os.path.dirname(__file__)), "PeleTemplates"))
        self.usesrun = "true" if args.usesrun else "false"

    def optative_params(self, args):
        self.input = args.input
        self.forcefield = args.forcefield if args.forcefield else self.simulation_params.get("forcefield", "OPLS2005")
        self.solvent = args.solvent if args.solvent else self.simulation_params.get("solvent", "VDGBNP") 
        self.verbose = "true" if args.verbose else self.simulation_params.get("verbose", "false")
        self.cpus = args.cpus = args.cpus if args.cpus else self.simulation_params.get("cpus", 60)
        self.restart = args.restart if args.restart else self.simulation_params.get("restart", "all")
        self.test = args.test 
        #+1 to avoid being 0 
        self.equil_steps = int(args.eq_steps/self.cpus) + 1 if args.eq_steps else self.simulation_params.get("equilibration_steps", 1)
        self.equilibration = "true" if args.equilibration else self.simulation_params.get("equilibration", "false")
        self.adaptive_restart = args.adaptive_restart
        self.poses = args.poses if args.poses else self.simulation_params.get("poses", self.cpus-1)
        self.pele_exec = args.pele_exec if args.pele_exec else os.path.join(cs.PELE, "bin/Pele_mpi")
        self.pele_data = args.pele_data if args.pele_data else os.path.join(cs.PELE, "Data")
        self.pele_documents = args.pele_documents if args.pele_documents else os.path.join(cs.PELE, "Documents")

    def system_preparation_params(self, args):
        self.skip_prep = args.skip_prep if args.skip_prep else self.simulation_params.get("skip_prep", False)
        self.nonstandard = args.nonstandard if args.nonstandard else self.simulation_params.get("nonstandard", [])
        self.constraints = None
        self.external_constraints = hp.retrieve_constraints_for_pele(args.external_constraints, self.system) if args.external_constraints else None
        self.constrain_smiles = args.constrain_smiles if args.constrain_smiles else self.simulation_params.get("constrain_smiles", None)
        self.no_ppp = args.no_ppp if args.no_ppp else self.simulation_params.get("no_ppp", False)

    def ligand_params(self, args):
        self.spython = os.path.join(cs.SCHRODINGER, "utilities/python")
        if not os.path.exists(self.spython):
            self.spython = os.path.join(cs.SCHRODINGER, "run")
        self.mae_lig = args.mae_lig if args.mae_lig else self.simulation_params.get("mae_lig", None)
        self.external_template = args.template if args.template else self.simulation_params.get("template", [])
        self.external_rotamers = args.rotamers if args.rotamers else self.simulation_params.get("rotamers", [])
        self.skip_ligand_prep = args.skip_ligand_prep if args.skip_ligand_prep else self.simulation_params.get("args.skip_ligand_prep", [])

    def water_params(self, args):
        self.water_temp = args.water_temp if args.water_temp else self.simulation_params.get("water_temp", 5000)
        self.water_overlap = args.water_overlap if args.water_overlap else self.simulation_params.get("water_overlap", 0.78)
        self.water_constr = args.water_constr if args.water_constr else self.simulation_params.get("water_constr", 0)
        self.water_trials = args.water_trials if args.water_trials  else self.simulation_params.get("water_trials", 10000)
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
        self.box_center = "["+ ",".join([str(coord) for coord in args.box_center]) + "]" if args.box_center else self.simulation_params.get("box_center", None)

    def metrics_params(self, args):
        self.metrics = ""
        self.native = cs.NATIVE.format(os.path.abspath(args.native), self.chain) if args.native else ""
        self.atom_dist = args.atom_dist
        self.metrics = mt.build_metrics(self.system, self.atom_dist)

    def output_params(self, args):
        self.folder = args.folder
        self.output = args.output if args.output is not None else "output"
        self.report_name = args.report_name if args.report_name else "report"
        self.traj_name = args.traj_name if args.traj_name else "trajectory.pdb"
        self.xtc = self.traj_name.endswith(".xtc")
        self.pdb = self.traj_name.endswith(".pdb")

    def analysis_params(self, args):
        self.analyse = args.analyse if args.analyse else True
        self.mae = args.mae if args.mae else False
        self.only_analysis = args.only_analysis
        self.analysis_nclust = args.analysis_nclust
        self.overwrite = args.overwrite
        self.be_column = args.be_column
        self.te_column = args.te_column
        self.limit_column = args.limit_column
