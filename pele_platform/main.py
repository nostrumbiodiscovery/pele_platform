import matplotlib
matplotlib.use("Agg")
import yaml
import sys
import pele_platform.constants.constants as cs
sys.path.append(cs.DIR)
from argparse import HelpFormatter
from operator import attrgetter
import argparse
import os
import pele_platform.Adaptive.simulation as ad
import pele_platform.Frag.simulation as fr
from pele_platform.PPI.main import run_ppi

class Launcher():


    def __init__(self, arguments):
        self.cpus = arguments.cpus
        self.restart = arguments.restart
        self.test = arguments.test
        self._args = arguments
        if arguments.frag_core:
            self.pele_feature = "frag"
        elif arguments.ppi:
            self.pele_feature = "PPI"
        else: 
            self.pele_feature = "adaptive"

    def launch(self):
        if not self._args.no_check:
            self._check_variables()
        if self.pele_feature == "adaptive":
            job_variables = ad.run_adaptive(self._args)
        elif self.pele_feature == "PPI":
            job_variables = run_ppi(self._args)
        elif self.pele_feature == "frag":
            #Set variables and input ready 
            job_variables = fr.FragRunner(self._args)
            job_variables.prepare_control_file()
            #Set test variables if desired
            if self.test:
                job_variables.set_test_variables()
            #Depending on input different method
            if job_variables.ligands: #Full ligands as sdf
                job_variables.prepare_input_file()
                job_variables.run()
            elif job_variables.ai:
                job_variables.grow_ai()
            else:
                job_variables.run()
            # Execute job
        return job_variables

    def _check_variables(self):
        variables = [ 
        EnvVariable("pele_data", self._args.pele_data, os.path.join(cs.PELE, "Data"), "--pele_data /path/to/data/folder/", "export PELE=/path/to/PELE-1.X/"),
        EnvVariable("pele_documents", self._args.pele_documents, os.path.join(cs.PELE, "Documents"), "--pele_documents /path/to/documents/folder", "export PELE=/path/to/PELE-1.X/"),
        EnvVariable("pele_exec", self._args.pele_exec, os.path.join(cs.PELE, "bin/Pele_mpi"), "--pele_exec /path/to/PELE_exec", "export PELE=/path/to/PELE-1.X/"),
        EnvVariable("pele_license", self._args.pele_license, os.path.join(cs.PELE, "licenses"), "--pele_license /path/to/licenses", "export PELE=/path/to/PELE-1.X/"),
        EnvVariable("schrodinger", self._args.schrodinger, cs.SCHRODINGER, "--schrodinger /path/to/schrodinger-20XX/", "export SCHRODINGER=/path/to/schrodinger-20XX/")
        ]
        for variable in variables:
            variable.is_valid() 


class EnvVariable():

    def __init__(self, name, variable, default, flag, env_var):
        self.name = name
        self.variable = variable
        self.default = default
        self.flag = flag
        self.env_var = env_var
        
    def is_valid(self):
        if self.variable:
            if os.path.exists(self.variable):
                return True
        elif self.default:
           if os.path.exists(self.default):
               return True

        raise ValueError("{} not found. If you have the standard installation \
export the environment variable by doing: {}.\
 Else define the location of the library via the flag: {}".format(self.name, self.env_var, self.flag))

class SortingHelpFormatter(HelpFormatter):
    def add_arguments(self, actions):
        actions = sorted(actions, key=attrgetter('option_strings'))
        super(SortingHelpFormatter, self).add_arguments(actions)


def parseargs_yaml(args=[]):
    parser = argparse.ArgumentParser(description='Run PELE Platform', formatter_class=SortingHelpFormatter)
    parser.add_argument('input_file', type=str, help='Yaml input file')
    args = parser.parse_args(args) if args else parser.parse_args()
    return args
    

class YamlParser(object):

    def __init__(self, yamlfile):
        self.yamlfile = yamlfile
        self.parse()

    def parse_yaml(self):
        with open(self.yamlfile, 'r') as stream:
            try:
                data = yaml.safe_load(stream)
            except yaml.YAMLError as exc:
                raise(exc)
        return data
    
    def parse(self):
        data = self.parse_yaml()
        self.system = data.get("system", "")
        self.system = os.path.abspath(self.system) if self.system else ""
        self.residue = data.get("resname", None)
        self.chain = data.get("chain", None)
        self.hbond = data.get("hbond", [None, None])
        self.test = data.get("test", None)
        self.pele = data.get("pele", None)
        self.forcefield = data.get("forcefield", "OPLS2005")
        self.verbose = data.get("verbose", None)
        self.anm_freq = data.get("anm_freq", None)
        self.sidechain_freq = data.get("sidechain_freq", None)
        self.min_freq = data.get("min_freq", None)
        self.water_freq = data.get("water_freq", None)
        self.temperature = self.temp = data.get("temperature", None)
        self.sidechain_resolution = data.get("sidechain_res", None)
        self.steric_trials = data.get("steric_trials", None)
        self.overlap_factor = data.get("overlap_factor", None)
        self.solvent = data.get("solvent", None)
        self.usesrun = data.get("usesrun", None)
        self.spawning = data.get("spawning", None)
        self.iterations = data.get("iterations", None)
        self.pele_steps = self.steps = data.get("steps", None)
        self.cpus = data.get("cpus", None)
        self.density = data.get("density", None)
        self.cluster_values = data.get("cluster_values", None)
        self.cluster_conditions = data.get("cluster_conditions", None)
        self.simulation_type = data.get("simulation_type", None)
        self.equilibration = data.get("equilibration", None)
        self.eq_steps = data.get("equilibration_steps", None)
        self.adaptive_restart = data.get("adaptive_restart", None)
        self.input = data.get("global_inputs", None)
        self.report_name = data.get("report", None)
        self.traj_name = data.get("traj", None)
        self.adaptive = data.get("adaptive", None)
        self.epsilon = data.get("epsilon", None)
        self.bias_column = data.get("bias_column", None)
        self.gridres = data.get("gridres", 10)
        self.core = data.get("core", -1)
        self.mtor = data.get("maxtorsion", 4)
        self.n = data.get("n", 10000)
        self.template = data.get("templates", None)
        self.ext_temp = self.template
        self.rotamers = data.get("rotamers", None)
        self.ext_rotamers = self.rotamers
        self.mae_lig = data.get("mae_lig", None)
        self.mae_lig = os.path.abspath(self.mae_lig) if self.mae_lig else None
        self.skip_prep = self.no_ppp = data.get("skip_preprocess", None)
        self.gaps_ter = data.get("TERs", None)
        self.charge_ter = data.get("charge_ters", None)
        self.nonstandard = data.get("nonstandard", None)
        self.prepwizard = data.get("prepwizard", None)
        self.box_center = data.get("box_center", None)
        self.box_center = [str(x) for x in self.box_center] if self.box_center else None
        self.box_radius = data.get("box_radius", None)
        self.box = data.get("box", None)
        self.native = data.get("rmsd_pdb", "")
        self.atom_dist = data.get("atom_dist", None)
        self.debug = data.get("debug", None)
        self.folder = data.get("working_folder", None)
        self.output = data.get("output", None)
        self.randomize = data.get("randomize", None)
        self.full = data.get("global", None)
        self.proximityDetection = data.get("proximityDetection", None)
        self.poses = data.get("poses", None)
        self.precision_glide = data.get("precision_glide", None) 
        self.msm = data.get("msm", None)
        self.precision = data.get("precision", None)
        self.clust = data.get("exit_clust", None)
        self.restart = data.get("msm_restart", None)
        self.lagtime = data.get("lagtime", None)
        self.msm_clust = data.get("msm_clust", None)
        self.rescoring = data.get("rescoring", None)
        self.in_out = data.get("in_out", None)
        self.in_out_soft = data.get("in_out_soft", None)
        self.exit = data.get("exit", None)
        self.exit_value = data.get("exit_value", None)
        self.exit_condition = data.get("exit_condition", None)
        self.exit_trajnum = data.get("exit_trajnum", None)
        self.water_exp = data.get("water_bs", None)
        self.water_lig = data.get("water_lig", None)
        self.water = data.get("water", None)
        self.water_expl = data.get("water_expl", None)
        self.water_freq = data.get("water_freq", None)
        self.water_center = data.get("box_water", None)
        self.water_temp = data.get("water_temp", None)
        self.water_overlap = data.get("water_overlap", None)
        self.water_constr = data.get("water_constr", None)
        self.water_trials = data.get("water_trials", None)
        self.water_radius = data.get("water_radius", None)
        self.bias = data.get("bias", None)
        self.induced_fit_exhaustive = data.get("induced_fit_exhaustive", None)
        self.induced_fit_fast = data.get("induced_fit_fast", None)
        self.frag = data.get("frag", None)
        self.ca_constr=data.get("ca_constr", None)
        self.one_exit=data.get("one_exit", None)
        self.box_type=data.get("box_type", None)
        self.box_metric = data.get("box_metric", None)
        self.time = data.get("time", None)
        self.nosasa = data.get("nosasa", None)
        self.sasa = data.get("sasa", None)
        self.perc_sasa = data.get("perc_sasa", None)
        self.seed=data.get("seed", None)
        self.pdb = data.get("pdb", None)
        self.log = data.get("log", None)
        self.nonrenum = data.get("nonrenum", None)
        self.pele_exec = data.get("pele_exec", None)
        self.pele_data = data.get("pele_data", None)
        self.pele_documents = data.get("pele_documents", None)
        self.pca = data.get("pca", None)
        self.anm_direction = data.get("anm_direction", None)
        self.anm_mix_modes = data.get("anm_mix_modes", None)
        self.anm_picking_mode = data.get("anm_picking_mode", None)
        self.anm_displacement = data.get("anm_displacement", None)
        self.anm_modes_change = data.get("anm_modes_change", None)
        self.anm_num_of_modes = data.get("anm_num_of_modes", None)
        self.anm_relaxation_constr = data.get("anm_relaxation_constr", None)
        self.remove_constraints = data.get("remove_constraints", None)
        self.pca_traj = data.get("pca_traj", None)
        self.perturbation = data.get("perturbation", None)
        self.binding_energy = data.get("binding_energy", None)
        self.sasa = data.get("sasa", None)
        self.parameters = data.get("parameters", None)
        self.analyse = data.get("analyse", None)
        self.selection_to_perturb = data.get("selection_to_perturb", None)
        self.mae = data.get("mae", None)
        self.constrain_smiles = data.get("constrain_smiles", None)
        self.skip_ligand_prep = data.get("skip_ligand_prep", None)
        self.spawning_condition = data.get("spawning_condition", None)
        self.external_constraints = data.get("external_constraints", None)
        self.only_analysis = data.get("only_analysis", False)
        self.overwrite = data.get("overwrite_analysis", True)
        self.analysis_nclust = data.get("analysis_nclust", 10)
        self.te_column = data.get("te_column", 4)
        self.be_column = data.get("be_column", 5)
        self.limit_column = data.get("limit_column", 6)
        self.com = data.get("COMligandConstraint", None)
        self.pele_license = data.get("pele_license", None)
        self.schrodinger = data.get("schrodinger", None)
        self.no_check = data.get("no_check", False)

        #Frag
        self.frag_core = data.get("frag_core", False)
        self.frag_input = data.get("frag_input", False)
        self.frag_ligands = data.get("frag_ligands", False)
        self.growing_steps = data.get("growing_steps", False)
        self.frag_steps = data.get("steps_in_gs", False)
        self.frag_eq_steps = data.get("sampling_steps", False)
        self.protocol = data.get("protocol", None)
        self.frag_ai = data.get("frag_ai", False)
        self.frag_ai_iterations = data.get("frag_ai_iterations", False)
        self.chain_core = data.get("chain_core", False)

        #PPI
        self.n_components = data.get("n_components", None)
        self.ppi = data.get("ppi", None)

        if self.test:
            print("##############################")
            print("WARNING: This simulation is a test do not use the input files to run production simulations")
            print("##############################")
            self.cpus = 5
            self.pele_steps = self.steps = 1
            self.iterations = 1
            self.min_freq = 0
            self.anm_freq = 0
            self.sidechain_freq = 0
            self.temperature = self.temp = 10000
            self.n_components = 3

def run_platform(input_yaml):
    arguments = parseargs_yaml([input_yaml,])
    arguments = YamlParser(arguments.input_file)
    job_params = Launcher(arguments).launch()
    return job_params

def run_platform(input_yaml):
    arguments = parseargs_yaml([input_yaml,])
    arguments = YamlParser(arguments.input_file)
    job_params = Launcher(arguments).launch()
    return job_params


if __name__ == "__main__":
    arguments = parseargs_yaml()
    arguments = YamlParser(arguments.input_file)
    job_params = Launcher(arguments).launch()
