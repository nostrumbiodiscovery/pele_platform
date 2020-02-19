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
import pele_platform.Adaptive.main as gl
import pele_platform.Adaptive.simulation as ad


class Launcher():


    def __init__(self, arguments):
        self.cpus = arguments.cpus
        self.restart = arguments.restart
        self.test = arguments.test
        self._args = arguments

    def launch(self):
        ad.run_adaptive(self._args)

class SortingHelpFormatter(HelpFormatter):
    def add_arguments(self, actions):
        actions = sorted(actions, key=attrgetter('option_strings'))
        super(SortingHelpFormatter, self).add_arguments(actions)

def parseargs(args=[]):
    parser = argparse.ArgumentParser(description='Run PELE Platform', formatter_class=SortingHelpFormatter)
    parser.add_argument('system', type=str, help='complex to run pele on')
    parser.add_argument('residue', type=str, help='residue of the ligand to extract', default=cs.LIG_RES)
    parser.add_argument('chain', type=str, help='chain of the ligand to extract', default=cs.LIG_CHAIN)
    parser.add_argument("--test", action='store_true', help="Run a fast pele_platform test")

    #Pele info
    pele_info = parser.add_argument_group('\nPele Parameters') 
    pele_info.add_argument("--pele", type=str,  help="Pele control_file", default=None)
    pele_info.add_argument('--forcefield', '-f', type=str, help='chain of the ligand to extract', default=cs.FORCEFIELD)
    pele_info.add_argument("--anm_freq", type=int,  help="Frequency to perform ANM", default=4)
    pele_info.add_argument("--sidechain_freq", type=int,  help="Frequency to perform sidechain sampling", default=2)
    pele_info.add_argument("--min_freq", type=int,  help="Frequency to perform all atoms minimization", default=1)
    pele_info.add_argument("--temperature", type=int,  help="Temperature to perform PELE simulation", default=1500)
    pele_info.add_argument("--sidechain_resolution", type=int,  help="Every how many degrees the sidechains will be rotated [10, 30...]", default=10)
    pele_info.add_argument("--steric_trials", type=int,  help="Number fo steric trials on sidechain sampling", default=None)
    pele_info.add_argument("--overlap_factor", type=float,  help="Relaxation of vanderwals clashes from 0 to 1", default=None)
    pele_info.add_argument('--solvent', type=str, help='Type of implicit solvent (OBC/VDGBNP). default [OBC]. i.e. --solvent VDGBNP', default="VDGBNP")

    # Adaptive info
    adaptive_info = parser.add_argument_group('\nAdaptive Parameters') 
    adaptive_info.add_argument("--usesrun", action="store_true",  help="Use srun option with intel")
    adaptive_info.add_argument("--iterations", type=int,  help="PELE iterations", default=50)
    adaptive_info.add_argument("--pele_steps", type=int,  help="PELE steps", default=None)
    adaptive_info.add_argument("--spawning", type=str,  help="Adaptive Spawning type. [epsilon, inverselyProportional]", default=None)
    adaptive_info.add_argument("--density", type=str,  help="Adaptive Density type. [null, continuous]", default=None)
    adaptive_info.add_argument("--cluster_values", type=str,  help="Cluster values", default=None)
    adaptive_info.add_argument("--cluster_conditions", type=str,  help="Cluster conditions", default=None)
    adaptive_info.add_argument("--simulation_type", type=str,  help="Type of simulation [pele, md]", default=None)
    adaptive_info.add_argument('--cpus', type=int, help='number of processors', default=cs.CPUS)
    adaptive_info.add_argument("--equilibration", action="store_true",  help="Whether to do a initial equilibration or not")
    adaptive_info.add_argument("--eq_steps", type=float,  help="Number of equilibration steps", default=cs.EQ_STEPS)
    adaptive_info.add_argument("--adaptive_restart", action="store_true",  help="Whether to set restart true to adaptive")
    adaptive_info.add_argument('--input', nargs="+", help='Set initial input for simulation')
    adaptive_info.add_argument('--report_name', type=str, help='Name of the output report files. i.e. report_Sim1_. Default [report]', default="report")
    adaptive_info.add_argument('--traj_name', type=str, help='Name of the trajectory files (pdb and xtc supported). i.e. trajectory_sim1.pdb. Default [trajectory.xtc]', default="trajectory.xtc")
    adaptive_info.add_argument("--adaptive", type=str,  help="Adaptive control_file", default=None)

    adaptive_info.add_argument("--epsilon", type=float,  help="From 0 to 1 how bias do you want the simulation", default=None)
    adaptive_info.add_argument("--bias_column", type=int,  help="Column of the report starting by 1 towards where to bias simulation", default=None)



    #Ligand Parametrization
    plop_info = parser.add_argument_group('\nLigand Parametrization') 
    plop_info.add_argument("--gridres", type=str, help="Rotamers angle resolution", default=cs.GRIDRES)
    plop_info.add_argument("--core", type=int, help="Give one atom of the core section", default=-1)
    plop_info.add_argument("--mtor", type=int, help="Gives the maximum number of torsions allowed in each group.  Will freeze bonds to extend the core if necessary.", default=4)
    plop_info.add_argument("--n", type=int, help="Maximum Number of Entries in Rotamer File", default=1000)
    plop_info.add_argument("--template", type=str,  help="External template for ligand", default=None)
    plop_info.add_argument("--rotamers", type=str,  help="External romtamers library for ligand", default=None)
    plop_info.add_argument("--mae_lig", type=str, help="ligand .mae file to include QM charges coming from jaguar")

    #PDB preparation info
    receptor_preparation = parser.add_argument_group('\nReceptor Preparation') 
    receptor_preparation.add_argument("--no_ppp", action="store_true",  help="Do not run experimental PPP")
    receptor_preparation.add_argument("--gaps_ter", help="Include TER when a possible gap is found", action='store_true')
    receptor_preparation.add_argument("--charge_ter", help="Charge protein terminals", action='store_true')
    receptor_preparation.add_argument("--skip_prep", action="store_true",  help="Whether to do the initial preprocessing or not")
    receptor_preparation.add_argument("--nonstandard", nargs="+",  help="Mid Chain non standard residues to be treated as ATOM not HETATOM", default = [])
    receptor_preparation.add_argument("--prepwizard", action="store_true",  help="Run protein preparation wizard from schrodinger")


    #Box info
    box_info = parser.add_argument_group('\nBox') 
    box_info.add_argument("--box", type=str, help="Exploration box for Pele")
    box_info.add_argument("--user_center", "-c", nargs='+', type=str, help='center of the box', default=None)
    box_info.add_argument("--box_radius", "-r", type=float,  help="Radius of the box", default=None)

    #Metrics info
    metrics_info = parser.add_argument_group('\nMetrics') 
    metrics_info.add_argument('--native', type=str, help='native file to compare RMSD to', default="")
    metrics_info.add_argument("--atom_dist", nargs="+",  help="Number of the atoms to calculate the distance in between i.e --atom dist 123 456 122 121", default=None)

    #Output info
    output_info = parser.add_argument_group('\nOutput info') 
    output_info.add_argument("--folder", "-wf", type=str,  help="Name of the output working folder", default=None)
    output_info.add_argument("--output", "-o", type=str,  help="Name of the adaptive output", default="output")

    #Full simulation
    full_info = parser.add_argument_group('\nFull Simulation') 
    full_info.add_argument('--randomize', action="store_true", help='Randomize ligand position around protein')
    full_info.add_argument("--full", action="store_true",  help="Launch ful ful binding site exploration")
    full_info.add_argument("--poses", type=int,  help="Number of ligand poses for global exploration", default=40)

    #Kinase simulation
    kinase_info = parser.add_argument_group('\nFull Simulation') 
    kinase_info.add_argument("--hbond", nargs='+',  help="Definition of kinase hbond", default="" )
    kinase_info.add_argument("--precision_glide", type=str,  help="Glide precision.. Options = [SP, XP]", default="SP")

    #MSM simulation
    msm_info = parser.add_argument_group('\nMSM Simulation') 
    msm_info.add_argument("--msm", action="store_true",  help="Launch MSM")
    msm_info.add_argument("--clust", type=int, help="Numbers of clusters to start PELE's exploration with", default=cs.CLUSTERS)
    msm_info.add_argument("--restart", type=str, help="Restart the platform from [all, pele, msm] with these keywords", default=cs.PLATFORM_RESTART)
    msm_info.add_argument("--precision", action='store_true', help="Use a more agressive control file to achieve better convergence")
    msm_info.add_argument("--lagtime", type=int,  help="MSM Lagtime", default=100)
    msm_info.add_argument("--msm_clust", type=int,  help="MSM cluster number", default=200)


    #in_out simulation
    dissociation_info = parser.add_argument_group('\nDissociation path Simulation') 
    dissociation_info.add_argument("--in_out", action="store_true",  help="Launch dissotiation path adaptive")
    dissociation_info.add_argument("--in_out_soft", action="store_true",  help="Launch soft dissotiation path adaptive")


    #Exit simulation
    exit_info = parser.add_argument_group('\nExit Simulation') 
    exit_info.add_argument("--exit", action="store_true",  help="Exit simulation given certain condition")
    exit_info.add_argument("--exit_value", type=float,  help="Value of the metric used to exit simulation", default=0.95)
    exit_info.add_argument("--exit_condition", type=str,  help="Selects wether to exit the simulation when being over o below the exit value", default=">")
    exit_info.add_argument("--exit_trajnum", type=str,  help="Number of trajectories to accomplished the condition to exit the simulation", default=4)

    #Water simulation
    water_info = parser.add_argument_group('\nWater Simulation') 
    water_info.add_argument("--water_exp", type=str,  help="Launch water exploration adaptive PELE", default=None)
    water_info.add_argument("--water_lig", nargs="+",  help="Launch ligand-water exploration adaptive PELE", default=None)
    water_info.add_argument("--water_center", nargs="+",  help="Launch ligand-water exploration adaptive PELE", default=None)
    water_info.add_argument("--water_temp",  type=int, help="Temperature of water MC. i.e --water_temp 1000 default=500", default=500)
    water_info.add_argument("--water_constr",  type=float, help="Constraint on the waters MC. i.e ----water_const 0.5 default=0.2", default=0.2)
    water_info.add_argument("--water_trials",  type=int, help="Steric trials on the waters MC. i.e --water_trials 2000 default=1000", default=1000)


    #Bias simulation
    bias_info = parser.add_argument_group('\nBias Simulation') 
    bias_info.add_argument("--bias", action="store_true",  help="Launch biased simulation")

    
    #Induced fit simulation
    induced_fit_info = parser.add_argument_group('\nInduced fit simulation') 
    induced_fit_info.add_argument("--induce_fit", action="store_true",  help="Launch induce fit adaptive")



    args = parser.parse_args(args) if args else parser.parse_args()

    return args


def parseargs_yaml(args=[]):
    parser = argparse.ArgumentParser(description='Run PELE Platform', formatter_class=SortingHelpFormatter)
    parser.add_argument('input_file', type=str, help='Yaml input file')
    args = parser.parse_args(args) if args else parser.parse_args()
    return args
    
def main(arguments):
    """
    Main function that sets the functionality
    of the software that will be used [Pele, Adaptive, glide...]
    and launch the respective job
    """
    job = Launcher(arguments)
    job.launch()
    return job



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
        self.system = os.path.abspath(data.get("system", None))
        self.residue = data.get("resname", None)
        self.chain = data.get("chain", None)
        self.hbond = data.get("hbond", [None, None])
        self.test = data.get("test", None)
        self.pele = data.get("pele", None)
        self.forcefield = data.get("forcefield", "OPLS2005")
        self.verbose = data.get("verbose", False)
        self.anm_freq = data.get("anm_freq", 4)
        self.sidechain_freq = data.get("sidechain_freq", 2)
        self.min_freq = data.get("min_freq", 1)
        self.water_freq = data.get("water_freq", 1)
        self.temperature = self.temp = data.get("temperature", 1500)
        self.sidechain_resolution = data.get("sidechain_res", 10)
        self.steric_trials = data.get("steric_trials", 250)
        self.overlap_factor = data.get("overlap_factor", 0.65)
        self.solvent = data.get("solvent", "VDGBNP")
        self.usesrun = data.get("usesrun", False)
        self.spawning = data.get("spawning", None)
        self.iterations = data.get("iterations", None)
        self.pele_steps = self.steps = data.get("steps", None)
        self.cpus = data.get("cpus", 2)
        self.density = data.get("density", None)
        self.cluster_values = data.get("cluster_values", None)
        self.cluster_conditions = data.get("cluster_conditions", None)
        self.simulation_type = data.get("simulation_type", None)
        self.equilibration = data.get("equilibration", False)
        self.eq_steps = data.get("equilibration_steps", 2.0)
        self.adaptive_restart = data.get("adaptive_restart", False)
        self.input = data.get("input", [])
        self.report_name = data.get("report", "report")
        self.traj_name = data.get("traj", "trajectory.pdb")
        self.adaptive = data.get("adaptive", None)
        self.epsilon = data.get("epsilon", None)
        self.bias_column = data.get("bias_column", None)
        self.gridres = data.get("gridres", 10)
        self.core = data.get("core", -1)
        self.mtor = data.get("maxtorsion", 4)
        self.n = data.get("n", 100000)
        self.template = data.get("templates", [])
        self.ext_temp = self.template
        self.rotamers = data.get("rotamers", [])
        self.ext_rotamers = self.rotamers
        self.mae_lig = data.get("mae_lig", None)
        self.mae_lig = os.path.abspath(self.mae_lig) if self.mae_lig else None
        self.skip_prep = self.no_ppp = data.get("preprocess", False)
        self.gaps_ter = data.get("TERs", False)
        self.charge_ter = data.get("charge_ters", False)
        self.nonstandard = data.get("nonstandard", [])
        self.prepwizard = data.get("prepwizard", False)
        self.user_center = data.get("box_center", None)
        self.user_center = [str(x) for x in self.user_center] if self.user_center else None
        self.box_radius = data.get("box_radius", None)
        self.box = data.get("box", None)
        self.native = data.get("rmsd_pdb", "")
        self.atom_dist = data.get("atom_dist", None)
        self.debug = data.get("debug", False)
        self.folder = data.get("working_folder", None)
        self.output = data.get("output", "output")
        self.randomize = data.get("randomize", False)
        self.full = data.get("global", False)
        self.proximityDetection = data.get("proximityDetection", True)
        self.poses = data.get("poses", 40)
        self.precision_glide = data.get("precision_glide", "SP") 
        self.msm = data.get("msm", False)
        self.precision = data.get("precision", False)
        self.clust = data.get("exit_clust", 40)
        self.restart = data.get("msm_restart", "all")
        self.lagtime = data.get("lagtime", 100)
        self.msm_clust = data.get("msm_clust", 200)
        self.rescoring = data.get("rescoring", False)
        self.in_out = data.get("in_out", False)
        self.in_out_soft = data.get("in_out_soft", False)
        self.exit = data.get("exit", False)
        self.exit_value = data.get("exit_value", 0.9)
        self.exit_condition = data.get("exit_condition", ">")
        self.exit_trajnum = data.get("exit_trajnum", 4)
        self.water_exp = data.get("water_bs", None)
        self.water_lig = data.get("water_lig", None)
        self.water = data.get("water", None)
        self.water_expl = data.get("water_expl", False)
        self.water_freq = data.get("water_freq", 1)
        self.water_center = data.get("box_water", None)
        self.water_temp = data.get("water_temp", 5000)
        self.water_overlap = data.get("water_overlap", 0.78)
        self.water_constr = data.get("water_constr", 0)
        self.water_trials = data.get("water_trials", 10000)
        self.water_radius = data.get("water_radius", None)
        self.bias = data.get("bias", False)
        self.induced_fit_exhaustive = data.get("induced_fit_exhaustive", False)
        self.induced_fit_fast = data.get("induced_fit_fast", False)
        self.frag = data.get("frag", False)
        self.ca_constr=data.get("ca_constr", 0.5)
        self.one_exit=data.get("one_exit", False)
        self.box_type=data.get("box_type", "fixed")
        self.box_metric = data.get("box_metric", None)
        self.time = data.get("time", False)
        self.nosasa = data.get("nosasa", False)
        self.sasa = data.get("sasa", False)
        self.perc_sasa = data.get("perc_sasa", [0.25, 0.5, 0.25])
        self.seed=data.get("seed", None)
        self.pdb = data.get("pdb", False)
        self.log = data.get("log", False)
        self.nonrenum = data.get("nonrenum", False)
        self.pele_exec = data.get("pele_exec", "")
        self.pele_data = data.get("pele_data", "")
        self.pele_documents = data.get("pele_documents", "")
        self.pca = data.get("pca", "")
        self.anm_direction = data.get("anm_direction", None)
        self.anm_mix_modes = data.get("anm_mix_modes", None)
        self.anm_picking_mode = data.get("anm_picking_mode", None)
        self.anm_displacement = data.get("anm_displacement", None)
        self.anm_modes_change = data.get("anm_modes_change", None)
        self.anm_num_of_modes = data.get("anm_num_of_modes", None)
        self.anm_relaxation_constr = data.get("anm_relaxation_constr", None)
        self.remove_constraints = data.get("remove_constraints", False)
        self.pca_traj = data.get("pca_traj", None)
        self.perturbation = data.get("perturbation", cs.PERTURBATION)
        self.binding_energy = data.get("binding_energy", cs.BE)
        self.sasa = data.get("sasa", cs.SASA)
        self.parameters = data.get("parameters", True)
        self.analyse = data.get("analyse", True)
        self.selection_to_perturb = data.get("selection_to_perturb", cs.SELECTION_TO_PERTURB)
        self.mae = data.get("mae", False)
        self.constrain_smiles = data.get("constrain_smiles", False)
        self.skip_ligand_prep = data.get("skip_ligand_prep", [])

        if self.test:
            print("##############################")
            print("WARNING: This simulation is a test do not use the input files to run production simulations")
            print("##############################")
            self.cpus = 2 if not self.full else 5
            self.pele_steps = self.steps = 1
            self.iterations = 1
            self.min_freq = 0
            self.anm_freq = 0
            self.sidechain_freq = 0
            self.temperature = self.temp = 10000


if __name__ == "__main__":
    arguments = parseargs_yaml()
    arguments = YamlParser(arguments.input_file)
    job = main(arguments)
