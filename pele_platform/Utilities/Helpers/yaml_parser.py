import os
import yaml




class YamlParser(object):

    def __init__(self, yamlfile):
        self.yamlfile = yamlfile
        self.checker()
        self.parse()

    def parse_yaml(self):
        with open(self.yamlfile, 'r') as stream:
            try:
                data = yaml.safe_load(stream)
            except yaml.YAMLError as exc:
                raise(exc)
        return data

    def checker(self):

        data = self.parse_yaml()

        valid_flags = {"system": "system",
        "residue": "resname",
        "chain": "chain",
        "hbond": "hbond",
        "test": "test",
        "pele": "pele",
        "forcefield": "forcefield",
        "verbose": "verbose",
        "anm_freq": "anm_freq",
        "sidechain_freq": "sidechain_freq",
        "min_freq": "min_freq",
        "water_freq": "water_freq",
        "temperature": "temperature",
        "sidechain_resolution": "sidechain_res",
        "steric_trials": "steric_trials",
        "overlap_factor": "overlap_factor",
        "solvent": "solvent",
        "usesrun": "usesrun",
        "spawning": "spawning",
        "iterations": "iterations",
        "pele_steps": "steps",
        "cpus": "cpus",
        "density": "density",
        "cluster_values": "cluster_values",
        "cluster_conditions": "cluster_conditions",
        "simulation_type": "simulation_type",
        "equilibration": "equilibration",
        "eq_steps": "equilibration_steps",
        "adaptive_restart": "adaptive_restart",
        "input": "global_inputs",
        "report_name": "report",
        "traj_name": "traj",
        "adaptive": "adaptive",
        "epsilon": "epsilon",
        "bias_column": "bias_column",
        "gridres": "gridres",
        "core": "core",
        "mtor": "maxtorsion",
        "n": "n",
        "template": "templates",
        "ext_temp": "template",
        "rotamers": "rotamers",
        "mae_lig": "mae_lig",
        "skip_prep": "skip_preprocess",
        "gaps_ter": "TERs",
        "charge_ter": "charge_ters",
        "nonstandard": "nonstandard",
        "prepwizard": "prepwizard",
        "box_center": "box_center",
        "box_radius": "box_radius",
        "box": "box",
        "native": "rmsd_pdb",
        "atom_dist": "atom_dist",
        "debug": "debug",
        "folder": "working_folder",
        "output": "output",
        "randomize": "randomize",
        "full": "global",
        "proximityDetection": "proximityDetection",
        "poses": "poses",
        "precision_glide": "precision_glide",
        "msm": "msm",
        "precision": "precision",
        "clust": "exit_clust",
        "restart": "restart",
        "lagtime": "lagtime",
        "msm_clust": "msm_clust",
        "rescoring": "rescoring",
        "in_out": "in_out",
        "in_out_soft": "in_out_soft",
        "exit": "exit",
        "exit_value": "exit_value",
        "exit_condition": "exit_condition",
        "exit_trajnum": "exit_trajnum",
        "water_exp": "water_bs",
        "water_lig": "water_lig",
        "water": "water",
        "water_expl": "water_expl",
        "water_freq": "water_freq",
        "water_center": "box_water",
        "water_temp": "water_temp",
        "water_overlap": "water_overlap",
        "water_constr": "water_constr",
        "water_trials": "water_trials",
        "water_radius": "water_radius",
        "bias": "bias",
        "induced_fit_exhaustive": "induced_fit_exhaustive",
        "induced_fit_fast": "induced_fit_fast",
        "frag": "frag",
        "ca_constr": "ca_constr",
        "one_exit": "one_exit",
        "box_type": "box_type",
        "box_metric": "box_metric",
        "time": "time",
        "nosasa": "nosasa",
        "sasa": "sasa",
        "perc_sasa": "perc_sasa",
        "seed": "seed",
        "pdb": "pdb",
        "log": "log",
        "nonrenum": "nonrenum",
        "pele_exec": "pele_exec",
        "pele_data": "pele_data",
        "pele_documents": "pele_documents",
        "pca": "pca",
        "anm_direction": "anm_direction",
        "anm_mix_modes": "anm_mix_modes",
        "anm_picking_mode": "anm_picking_mode",
        "anm_displacement": "anm_displacement",
        "anm_modes_change": "anm_modes_change",
        "anm_num_of_modes": "anm_num_of_modes",
        "anm_relaxation_constr": "anm_relaxation_constr",
        "remove_constraints": "remove_constraints",
        "pca_traj": "pca_traj",
        "perturbation": "perturbation",
        "binding_energy": "binding_energy",
        "sasa": "sasa",
        "parameters": "parameters",
        "analyse": "analyse",
        "selection_to_perturb": "selection_to_perturb",
        "mae": "mae",
        "constrain_smiles": "constrain_smiles",
        "skip_ligand_prep": "skip_ligand_prep",
        "spawning_condition": "spawning_condition",
        "external_constraints": "external_constraints",
        "only_analysis": "only_analysis",
        "overwrite": "overwrite_analysis",
        "analysis_nclust": "analysis_nclust",
        "te_column": "te_column",
        "be_column": "be_column",
        "limit_column": "limit_column",
        "com": "COMligandConstraint",
        "pele_license": "pele_license",
        "schrodinger": "schrodinger",
        "no_check": "no_check",
        "frag_core": "frag_core",
        "frag_input": "frag_input",
        "frag_ligands": "frag_ligands",
        "growing_steps": "growing_steps",
        "frag_steps": "steps_in_gs",
        "frag_eq_steps": "sampling_steps",
        "protocol": "protocol",
        "frag_ai": "frag_ai",
        "frag_ai_iterations": "frag_ai_iterations",
        "frag_run": "frag_run",
        "frag_restart": "frag_restart",
        "frag_output_folder": "frag_output_folder", 
        "chain_core": "chain_core",
        "n_components": "n_components",
        "frag_criteria": "frag_criteria",
        "frag_cluster_folder": "frag_cluster_folder",
        "ppi": "ppi",
        "rna": "rna"}

        for key in data.keys():
            if key not in valid_flags.values():
                raise KeyError("Input file contains an invalid keyword: {}".format(key))
        
        return valid_flags     
    
    def parse(self):
        valid_flags = self.checker()
        data = self.parse_yaml()
        self.system = data.get(valid_flags["system"], "")
        self.system = os.path.abspath(self.system) if self.system else ""
        self.residue = data.get(valid_flags["residue"], None)
        self.chain = data.get(valid_flags["chain"], None)
        self.hbond = data.get(valid_flags["hbond"], [None, None])
        self.test = data.get(valid_flags["test"], None)
        self.pele = data.get(valid_flags["pele"], None)
        self.forcefield = data.get(valid_flags["forcefield"], "OPLS2005")
        self.verbose = data.get(valid_flags["verbose"], None)
        self.anm_freq = data.get(valid_flags["anm_freq"], None)
        self.sidechain_freq = data.get(valid_flags["sidechain_freq"], None)
        self.min_freq = data.get(valid_flags["min_freq"], None)
        self.water_freq = data.get(valid_flags["water_freq"], None)
        self.temperature = self.temp = data.get(valid_flags["temperature"], None)
        self.sidechain_resolution = data.get(valid_flags["sidechain_resolution"], None)
        self.steric_trials = data.get(valid_flags["steric_trials"], None)
        self.overlap_factor = data.get(valid_flags["overlap_factor"], None)
        self.solvent = data.get(valid_flags["solvent"], None)
        self.usesrun = data.get(valid_flags["usesrun"], None)
        self.spawning = data.get(valid_flags["spawning"], None)
        self.iterations = data.get(valid_flags["iterations"], None)
        self.pele_steps = self.steps = data.get(valid_flags["pele_steps"], None)
        self.cpus = data.get(valid_flags["cpus"], None)
        self.density = data.get(valid_flags["density"], None)
        self.cluster_values = data.get(valid_flags["cluster_values"], None)
        self.cluster_conditions = data.get(valid_flags["cluster_conditions"], None)
        self.simulation_type = data.get(valid_flags["simulation_type"], None)
        self.equilibration = data.get(valid_flags["equilibration"], None)
        self.eq_steps = data.get(valid_flags["eq_steps"], None)
        self.adaptive_restart = data.get(valid_flags["adaptive_restart"], None)
        self.input = data.get(valid_flags["input"], None)
        self.report_name = data.get(valid_flags["report_name"], None)
        self.traj_name = data.get(valid_flags["traj_name"], None)
        self.adaptive = data.get(valid_flags["adaptive"], None)
        self.epsilon = data.get(valid_flags["epsilon"], None)
        self.bias_column = data.get(valid_flags["bias_column"], None)
        self.gridres = data.get(valid_flags["gridres"], 10)
        self.core = data.get(valid_flags["core"], -1)
        self.mtor = data.get(valid_flags["mtor"], 4)
        self.n = data.get(valid_flags["n"], 10000)
        self.template = data.get(valid_flags["template"], None)
        self.ext_temp = self.template
        self.rotamers = data.get(valid_flags["rotamers"], None)
        self.ext_rotamers = self.rotamers
        self.mae_lig = data.get(valid_flags["mae_lig"], None)
        self.mae_lig = os.path.abspath(self.mae_lig) if self.mae_lig else None
        self.skip_prep = self.no_ppp = data.get(valid_flags["skip_prep"], None)
        self.gaps_ter = data.get(valid_flags["gaps_ter"], None)
        self.charge_ter = data.get(valid_flags["charge_ter"], None)
        self.nonstandard = data.get(valid_flags["nonstandard"], None)
        self.prepwizard = data.get(valid_flags["prepwizard"], None)
        self.box_center = data.get(valid_flags["box_center"], None)
        self.box_center = [str(x) for x in self.box_center] if self.box_center else None
        self.box_radius = data.get(valid_flags["box_radius"], None)
        self.box = data.get(valid_flags["box"], None)
        self.native = data.get(valid_flags["native"], "")
        self.atom_dist = data.get(valid_flags["atom_dist"], None)
        self.debug = data.get(valid_flags["debug"], None)
        self.folder = data.get(valid_flags["folder"], None)
        self.output = data.get(valid_flags["output"], None)
        self.randomize = data.get(valid_flags["randomize"], None)
        self.full = data.get(valid_flags["full"], None)
        self.proximityDetection = data.get(valid_flags["proximityDetection"], None)
        self.poses = data.get(valid_flags["poses"], None)
        self.precision_glide = data.get(valid_flags["precision_glide"], None) 
        self.msm = data.get(valid_flags["msm"], None)
        self.precision = data.get(valid_flags["precision"], None)
        self.clust = data.get(valid_flags["clust"], None)
        self.restart = data.get(valid_flags["restart"], None)
        self.lagtime = data.get(valid_flags["lagtime"], None)
        self.msm_clust = data.get(valid_flags["msm_clust"], None)
        self.rescoring = data.get(valid_flags["rescoring"], None)
        self.in_out = data.get(valid_flags["in_out"], None)
        self.in_out_soft = data.get(valid_flags["in_out_soft"], None)
        self.exit = data.get(valid_flags["exit"], None)
        self.exit_value = data.get(valid_flags["exit_value"], None)
        self.exit_condition = data.get(valid_flags["exit_condition"], None)
        self.exit_trajnum = data.get(valid_flags["exit_trajnum"], None)
        self.water_exp = data.get(valid_flags["water_exp"], None)
        self.water_lig = data.get(valid_flags["water_lig"], None)
        self.water = data.get(valid_flags["water"], None)
        self.water_expl = data.get(valid_flags["water_expl"], None)
        self.water_freq = data.get(valid_flags["water_freq"], None)
        self.water_center = data.get(valid_flags["water_center"], None)
        self.water_temp = data.get(valid_flags["water_temp"], None)
        self.water_overlap = data.get(valid_flags["water_overlap"], None)
        self.water_constr = data.get(valid_flags["water_constr"], None)
        self.water_trials = data.get(valid_flags["water_trials"], None)
        self.water_radius = data.get(valid_flags["water_radius"], None)
        self.bias = data.get(valid_flags["bias"], None)
        self.induced_fit_exhaustive = data.get(valid_flags["induced_fit_exhaustive"], None)
        self.induced_fit_fast = data.get(valid_flags["induced_fit_fast"], None)
        self.frag = data.get(valid_flags["frag"], None)
        self.ca_constr=data.get(valid_flags["ca_constr"], None)
        self.one_exit=data.get(valid_flags["one_exit"], None)
        self.box_type=data.get(valid_flags["box_type"], None)
        self.box_metric = data.get(valid_flags["box_metric"], None)
        self.time = data.get(valid_flags["time"], None)
        self.nosasa = data.get(valid_flags["nosasa"], None)
        self.sasa = data.get(valid_flags["sasa"], None)
        self.perc_sasa = data.get(valid_flags["perc_sasa"], None)
        self.seed=data.get(valid_flags["seed"], None)
        self.pdb = data.get(valid_flags["pdb"], None)
        self.log = data.get(valid_flags["log"], None)
        self.nonrenum = data.get(valid_flags["nonrenum"], None)
        self.pele_exec = data.get(valid_flags["pele_exec"], None)
        self.pele_data = data.get(valid_flags["pele_data"], None)
        self.pele_documents = data.get(valid_flags["pele_documents"], None)
        self.pca = data.get(valid_flags["pca"], None)
        self.anm_direction = data.get(valid_flags["anm_direction"], None)
        self.anm_mix_modes = data.get(valid_flags["anm_mix_modes"], None)
        self.anm_picking_mode = data.get(valid_flags["anm_picking_mode"], None)
        self.anm_displacement = data.get(valid_flags["anm_displacement"], None)
        self.anm_modes_change = data.get(valid_flags["anm_modes_change"], None)
        self.anm_num_of_modes = data.get(valid_flags["anm_num_of_modes"], None)
        self.anm_relaxation_constr = data.get(valid_flags["anm_relaxation_constr"], None)
        self.remove_constraints = data.get(valid_flags["remove_constraints"], None)
        self.pca_traj = data.get(valid_flags["pca_traj"], None)
        self.perturbation = data.get(valid_flags["perturbation"], None)
        self.binding_energy = data.get(valid_flags["binding_energy"], None)
        self.parameters = data.get(valid_flags["parameters"], None)
        self.analyse = data.get(valid_flags["analyse"], None)
        self.selection_to_perturb = data.get(valid_flags["selection_to_perturb"], None)
        self.mae = data.get(valid_flags["mae"], None)
        self.constrain_smiles = data.get(valid_flags["constrain_smiles"], None)
        self.skip_ligand_prep = data.get(valid_flags["skip_ligand_prep"], None)
        self.spawning_condition = data.get(valid_flags["spawning_condition"], None)
        self.external_constraints = data.get(valid_flags["external_constraints"], None)
        self.only_analysis = data.get(valid_flags["only_analysis"], False)
        self.overwrite = data.get(valid_flags["overwrite"], True)
        self.analysis_nclust = data.get(valid_flags["analysis_nclust"], 10)
        self.te_column = data.get(valid_flags["te_column"], 4)
        self.be_column = data.get(valid_flags["be_column"], 5)
        self.limit_column = data.get(valid_flags["limit_column"], 6)
        self.com = data.get(valid_flags["com"], None)
        self.pele_license = data.get(valid_flags["pele_license"], None)
        self.schrodinger = data.get(valid_flags["schrodinger"], None)
        self.no_check = data.get(valid_flags["no_check"], False)

        #Frag
        self.frag_run = data.get(valid_flags["frag_run"], True)
        self.frag_core = data.get(valid_flags["frag_core"], False)
        self.frag_input = data.get(valid_flags["frag_input"], False)
        self.frag_ligands = data.get(valid_flags["frag_ligands"], False)
        self.growing_steps = data.get(valid_flags["growing_steps"], False)
        self.frag_steps = data.get(valid_flags["frag_steps"], False)
        self.frag_eq_steps = data.get(valid_flags["frag_eq_steps"], False)
        self.protocol = data.get(valid_flags["protocol"], None)
        self.frag_ai = data.get(valid_flags["frag_ai"], False)
        self.frag_ai_iterations = data.get(valid_flags["frag_ai_iterations"], False)
        self.chain_core = data.get(valid_flags["chain_core"], False)
        self.frag_restart = data.get(valid_flags["frag_restart"], False)
        self.frag_criteria = data.get(valid_flags["frag_criteria"], False)
        self.frag_output_folder = data.get(valid_flags["frag_output_folder"], False)
        self.frag_cluster_folder = data.get(valid_flags["frag_cluster_folder"], False)

        #PPI
        self.n_components = data.get(valid_flags["n_components"], None)
        self.ppi = data.get(valid_flags["ppi"], None)

        self.rna = data.get(valid_flags["rna"], None)

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

