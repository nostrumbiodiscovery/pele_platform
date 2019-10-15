import os
import shutil
import warnings
import random
import logging
import msm_pele.constants as cs
import msm_pele.Helpers.helpers as hp


class EnviroBuilder(object):
    """
        Base class wher the needed pele selfironment
        is build by creating folders and files
    """

    def __init__(self, folders, files, args):
        """
        Base class that encodes as attributes
        the software parameters for each stage
        """
        #Main Arguments
        self.folders = folders
        self.ext_temp = args.ext_temp
        self.files = files
        self.system = args.system
        self.box = args.box
        self.one_exit = args.one_exit
        self.user_center = args.user_center
        self.solvent = args.solvent
        self.user_radius = args.box_radius
        self.box_type = args.box_type
        self.box_metric = cs.BOX_METRIC if args.box_metric else " "
        self.iterations = args.iterations
        self.forcefield = args.forcefield
        self.water = args.water if args.water else None
        self.residue = args.residue
        self.templates = os.path.abspath(os.path.join(os.path.dirname(os.path.dirname(__file__)), "PeleTemplates"))
        self.restart = args.restart
        self.native = args.native
        self.chain = args.chain
        self.mae_lig = os.path.abspath(args.mae_lig) if args.mae_lig else None
        self.clusters = args.clust = args.clust if not args.test else 2
        self.exit_iters = args.exit_iters if not args.test else 1
        self.eq_struct = args.eq_struct if not args.test else 1
        self.test = args.test
        self.folder = args.folder
        self.pdb = args.pdb
        self.steps = args.steps
        self.nonstandard = args.nonstandard
        self.time = args.time
        self.lagtime = args.lagtime
        self.msm_clust = args.msm_clust
        self.log = '"simulationLogPath" : "$OUTPUT_PATH/logFile.txt",' if args.log else ""
        self.renumber = args.nonrenum
        self.nosasa = args.nosasa
        self.temp = args.temp
        self.sasa = args.sasa
        self.perc_sasa = args.perc_sasa
        self.water_radius = args.water_radius
        self.water_temp = args.water_temp
        self.water_trials = args.water_trials
        self.water_constr = args.water_constr
        self.usesrun = args.usesrun
        #Choose CPUS
        if args.test:
            self.cpus = args.cpus = 4
        elif args.restart == "analise":
            self.cpus = args.cpus = 1
        else:
            self.cpus = args.cpus
        #Build constants for each module
        self.build_sasa_constants()
        self.build_msm_constants()
        self.build_water_constants()
        self.build_path_constants()


    @classmethod
    def build_env(cls, args):
        if args.test and not args.precision:
            self = cls(cs.FOLDERS, cs.FILES_TEST, args)
        elif args.test and args.precision:
            self = cls(cs.FOLDERS, cs.FILES_TEST_XP, args)
        elif args.precision:
            self = cls(cs.FOLDERS, cs.FILES_XP, args)
        elif args.precision2:
            self = cls(cs.FOLDERS, cs.FILES_XP2, args)
        else:
            self = cls(cs.FOLDERS, cs.FILES_SP, args)
        self.create()
        return self


    def build_msm_constants(self):
        """
        Build sasa related constants for later
        classifing the exit simulation clusters
        """
        self.lagtime = 1 if self.test else self.lagtime
        self.lagtimes = None if self.test else [50, 100, 200, 500]
        self.msm_clust = 2 if self.test else self.msm_clust
        if not self.test:
            if self.time: self.steps = 10000
            else: self.steps = self.steps
        else:
            self.steps = 10


    def build_water_constants(self):
        """
        Build water constants to run
        water MC PELE
        """
        if self.water:
            cms = [ hp.find_coords(self.system, water.split(":")[1], water.split(":")[0]) for water in self.water]
            try:
                cm = [str(coord) for coord in hp.find_centroid(cms)]
            except TypeError:
                raise TypeError("Check the specified waters exist")
            water_atoms = [ '"' + water + '"' for water in self.water] 
            self.dynamic_water = cs.WATER.format(self.water_radius, ",".join(cm), ",".join(water_atoms),
                self.water_temp, self.water_trials, self.water_constr)
        else:
            self.water = []
            self.dynamic_water = ""

    def build_sasa_constants(self):
        """
        Build sasa related constants for later
        classifing the exit simulation clusters
        """
        self.perc_sasa_min, self.perc_sasa_int, self.perc_sasa_max = self.perc_sasa
        self.sasamin, self.sasamax = self.sasa if self.sasa else [None, None]
        self.sasa = True if not self.nosasa and not self.test else False


    def build_path_constants(self):

        self.template = None
        self.rotamers_file = None
        self.random_num = random.randrange(1, 70000)
        self.license = '''"{}"'''.format(cs.LICENSE)

        if self.test:
            self.equil_steps = 1
        else:
            self.equil_steps = int(cs.EQ_STEPS/self.cpus) if self.cpus < cs.EQ_STEPS else 1

        pele_dir = os.path.abspath("{}_Pele".format(self.residue))

        if not self.folder:
            self.pele_dir = is_repited(pele_dir) if self.restart == "all" else is_last(pele_dir)
        else:
            self.pele_dir = os.path.abspath(self.folder)

        if self.mae_lig:
            self.system_fix = os.path.join(self.pele_dir, "{}_complex_processed.pdb".format(os.path.splitext(os.path.basename(self.system))[0]))
        else:
            self.system_fix = os.path.join(self.pele_dir, "{}_processed.pdb".format(os.path.splitext(os.path.basename(self.system))[0]))

        for f in self.ext_temp:
            cs.FILES_NAME.append(os.path.join("DataLocal/Templates/{}/HeteroAtoms/".format(self.forcefield), os.path.basename(f)))
            self.files.append(os.path.basename(f))
            
        self.adap_ex_input = os.path.join(self.pele_dir, os.path.basename(self.system_fix))
        self.adap_exit_template = os.path.join(cs.DIR, "Templates/adaptive_exit.conf")
        self.adap_ex_output = os.path.join(self.pele_dir, "output_adaptive_exit")
        self.exit_path = os.path.join(self.adap_ex_output, "exit_path{}")
        self.template_folder = os.path.join(self.pele_dir, "DataLocal/Templates/{}/HeteroAtoms/".format(self.forcefield))
        self.obc_tmp = os.path.join(cs.DIR, "Templates/solventParamsHCTOBC.txt")
        self.obc_file = os.path.join(self.pele_dir, "DataLocal/OBC/solventParamsHCTOBC.txt")
        self.results = os.path.join(self.pele_dir, "results")
        self.cluster_output = os.path.join(self.pele_dir, "output_clustering")
        self.adap_l_input = os.path.join(self.pele_dir, "output_clustering/initial_*")
        self.adap_l_output = os.path.join(self.pele_dir, "output_pele")
        self.ad_ex_temp = os.path.join(self.pele_dir, "adaptive_exit.conf")
        self.ad_l_temp = os.path.join(self.pele_dir, "adaptive_long.conf")
        self.pele_exit_temp = os.path.join(self.pele_dir, "pele_exit.conf")
        self.pele_temp = os.path.join(self.pele_dir, "pele.conf")
        self.box_temp = os.path.join(self.pele_dir, "box.pdb")
        self._pele_temp = os.path.join(cs.DIR, "Templates/pele_SP.conf")
        self.clusters_output = os.path.join(self.cluster_output, "clusters_{}_KMeans_allSnapshots.pdb".format(self.clusters))
        self.ligand_ref = os.path.join(self.pele_dir, "ligand.pdb")
        if self.native:
            self.native = cs.NATIVE.format(os.path.abspath(self.native), self.chain)

        else:
            self.native = cs.NATIVE.format(os.path.abspath(self.ligand_ref), self.chain)
        self.topology = None if self.pdb else os.path.join(self.adap_ex_output, "topologies/topology_0.pdb")

    def update_variable_for_iteration(self, i):
        self.adap_ex_output = os.path.join(self.pele_dir, "output_adaptive_exit/iteration{}".format(i+1)) 
        self.topology = None if self.pdb else os.path.join(self.adap_ex_output, "topologies/topology_0.pdb")
        self.cluster_output = os.path.join(self.pele_dir, "output_clustering/iteration{}".format(i+1))
        self.clusters_output = os.path.join(self.cluster_output, "clusters_{}_KMeans_allSnapshots.pdb".format(self.clusters))
        self.adap_l_input = os.path.join(self.cluster_output, "initial_*")

    def create(self):
        if self.restart == "all":
            self.create_folders()
            self.create_files()
            self.create_logger()
        else:
            self.create_logger()

    def create_folders(self):
        """
            Create pele folders
        """

        for folder in self.folders:
            self.create_dir(self.pele_dir, folder)

    def create_files(self):
        """
            Copy templates
        """

        # Actions
        for file, destination_name in zip(self.files, cs.FILES_NAME):
            self.copy(file, os.path.join(self.pele_dir, destination_name))

    def create_dir(self, base_dir, extension=None):
        """
            Class Method to manage
            directory creation only if that
            ones doesn't exist

            Location:
                base_dir+extension
                or base_dir if extension is None
        """
        if extension:
            path = os.path.join(base_dir, extension)
            if os.path.isdir(path):
                warnings.warn("Directory {} already exists.".format(path), RuntimeWarning)
            else:
                os.makedirs(path)
        else:
            if os.path.isdir(base_dir):
                warnings.warn("Directory {} already exists.".format(base_dir), RuntimeWarning)
            else:
                os.makedirs(base_dir)

    def copy(self, standard, destination, user=None):
        if user:
            shutil.copy(user, os.path.join(self.pele_dir, standard))
        else:
            shutil.copy(standard, destination)
        return os.path.join(self.pele_dir, standard)

    def create_logger(self):
        log_name = os.path.join(self.pele_dir, "{}.log".format(self.residue))
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)
        formatter = logging.Formatter("%(asctime)s:%(levelname)s:%(message)s")
        if self.restart == "all":
            file_handler = logging.FileHandler(log_name, mode='w')
        else:
            file_handler = logging.FileHandler(log_name, mode='a')
        file_handler.setFormatter(formatter)
        self.logger.addHandler(file_handler)



def is_repited(pele_dir):

    original_dir = None
    split_dir = pele_dir.split("_")
    for chunk in split_dir:
        if chunk != "Pele":
            if original_dir:
                original_dir = "{}_{}".format(original_dir, chunk)
            else:
                original_dir = chunk
        else:
            break
    if split_dir[-1].isdigit():
        i = split_dir[-1]
        i = int(i) + 1 
    else:
        i = 1
    if os.path.isdir(pele_dir):
        new_pele_dir = "{}_Pele_{}".format(original_dir, i)
        new_pele_dir = is_repited(new_pele_dir)
        return new_pele_dir
    else:
        return pele_dir

def is_last(pele_dir):

    original_dir = None
    split_dir = pele_dir.split("_")
    for chunk in split_dir:
        if chunk != "Pele":
            if original_dir:
                original_dir = "{}_{}".format(original_dir, chunk)
            else:
                original_dir = chunk
        else:
            break
    if split_dir[-1].isdigit():
        i = split_dir[-1]
        i = int(i) + 1 
    else:
        i = 1 

    if os.path.isdir(pele_dir):
        new_pele_dir = "{}_Pele_{}".format(original_dir, i)
        if not os.path.isdir(new_pele_dir):
            return pele_dir
        else:
            new_pele_dir = is_last(new_pele_dir)
            return new_pele_dir
    else:
        return pele_dir

