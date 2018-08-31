import os
import shutil
import warnings
import random
import logging
import PELE_Platform.constants as cs


class EnviroBuilder(object):
    """
        Base class wher the needed pele environment
        is build by creating folders and files
    """

    def __init__(self, folders, files, file_names, args):
        self.folders = folders
        self.files = files
        self.file_names = file_names
        self.build_variables(args)
        self.build_paths()

    def build_variables(self, args):

        self.system = args.system
        self.forcefield = args.forcefield
        self.residue = args.residue
        self.templates = os.path.abspath(os.path.join(os.path.dirname(os.path.dirname(__file__)), "PeleTemplates"))
        self.cpus = args.cpus = args.cpus if not args.test else 4 
        self.restart = args.restart
        self.native = args.native
        self.chain = args.chain
        self.mae_lig = args.mae_lig
        self.clusters = args.clust = args.clust if not args.test else 2
        self.test = args.test
        self.folder = args.folder
        self.pdb = args.pdb
        self.hbond_donor, self.hbond_acceptor = args.hbond
        self.precision_glide = args.precision_glide
        self.adaptive = args.adaptive if args.adaptive else None
        self.pele = args.pele if args.adaptive else None
        self.software = args.software
        self.external_template = args.template
        self.external_rotamers = args.rotamers
        self.lagtime = args.lagtime
        self.msm_clust = args.msm_clust
        self.seed = random.randrange(1, 70000)
        self.license = '''"{}"'''.format(cs.LICENSE)
        self.equil_steps = 1 if self.test else int(cs.EQ_STEPS/self.cpus) + 1 #+1 to avoid being 0
        self.sasa_max = None
        self.sasa_min = None
        self.box_radius = None
        self.box_center = None
        self.constraints = None

    @classmethod
    def build_env(cls, args):
        FOLDERS, FILES, FILE_NAMES = cls.choose_env(args)
        env = cls(FOLDERS, FILES, FILE_NAMES, args)
        env.create()
        return env

    @staticmethod
    def choose_env(args):
        if args.test and args.software == "msm":
            return  cs.FOLDERS, cs.FILES_TEST, cs.FILES_NAME_MSM
        elif args.test and args.software == "glide":
            return cs.FOLDERS_GLIDE, cs.FILES_GLIDE_TEST, cs.FILES_NAME_GLIDE
        elif args.software == "glide" and not args.test:
            return cs.FOLDERS_GLIDE, cs.FILES_GLIDE, cs.FILES_NAME_GLIDE
        elif args.software == "msm" and args.precision:
            return  cs.FOLDERS, cs.FILES_XP, cs.FILES_NAME_MSM
        elif args.software == "msm" and not args.precision:
            return cs.FOLDERS, cs.FILES_SP, cs.FILES_NAME_MSM
        elif args.software == "adaptive":
            return cs.FOLDERS_ADAPTIVE, [args.adaptive, args.pele], [os.path.basename(args.adaptive), os.path.basename(args.pele)]
        elif args.software == "out_in":
            return cs.FOLDERS_ADAPTIVE, cs.FILES_OUT_IN, cs.FILES_NAME_OUT_IN
        elif args.software == "induce_fit":
            return cs.FOLDERS_ADAPTIVE, cs.FILES_INDUCE_FIT, cs.FILES_NAME_INDUCE_FIT
            
            
            
    def build_paths(self):


        pele_dir = os.path.abspath("{}_Pele".format(self.residue))

        if not self.folder:
            self.pele_dir = is_repited(pele_dir) if self.restart in cs.FIRST_RESTART else is_last(pele_dir)
        else:
            self.pele_dir = os.path.abspath(self.folder)

        if self.mae_lig:
            self.system_fix = os.path.join(self.pele_dir, "{}_complex_processed.pdb".format(os.path.abspath(os.path.splitext(self.system)[0])))
        else:
            self.system_fix = os.path.join(self.pele_dir, "{}_processed.pdb".format(os.path.abspath(os.path.splitext(self.system)[0])))

        self.rotamers_folder = os.path.join(self.pele_dir, "DataLocal/LigandRotamerLibs/")
        self.template_folder = os.path.join(self.pele_dir, "DataLocal/Templates/{}/HeteroAtoms/".format(self.forcefield))


        if self.software == "msm":
            self.adap_ex_input = os.path.join(self.pele_dir, os.path.basename(self.system_fix))
            self.adap_ex_output = os.path.join(self.pele_dir, "output_adaptive_exit")
            self.cluster_output = os.path.join(self.pele_dir, "output_clustering")
            self.adap_l_input = "{}/initial_*"
            self.adap_l_output = os.path.join(self.pele_dir, "output_pele")
            self.ad_ex_temp = os.path.join(self.pele_dir, "adaptive_exit.conf")
            self.ad_l_temp = os.path.join(self.pele_dir, "adaptive_long.conf")
            self.pele_exit_temp = os.path.join(self.pele_dir, "pele_exit.conf")
            self.pele_temp = os.path.join(self.pele_dir, "pele.conf")
            self.box_temp = os.path.join(self.pele_dir, "box.pdb")
            self.clusters_output = os.path.join(self.cluster_output, "clusters_{}_KMeans_allSnapshots.pdb".format(self.clusters))
            self.ligand_ref = os.path.join(self.pele_dir, "ligand.pdb")
            self.native = cs.NATIVE.format(os.path.abspath(self.native), self.chain) if self.native else cs.NATIVE.format(os.path.abspath(self.ligand_ref), self.chain)
            self.topology = None if self.pdb else os.path.join(self.adap_ex_output, "topology.pdb")

        elif self.software == "glide":
             self.adap_ex_input = os.path.join(self.pele_dir, os.path.basename(self.system_fix))
             self.adap_ex_output = os.path.join(self.pele_dir, "output_adaptive")
             self.cluster_output = os.path.join(self.pele_dir, "output_clustering")
             self.glide_template = os.path.join(self.pele_dir, "glide.in") 
             self.ad_ex_temp = os.path.join(self.pele_dir, "adaptive.conf")
             self.pele_exit_temp = os.path.join(self.pele_dir, "pele.conf")
             self.clusters_output = os.path.join(self.cluster_output, "clusters_{}_KMeans_allSnapshots.pdb".format(self.clusters))
             self.pele_temp = os.path.join(self.pele_dir, "pele.conf")
             self.ligand_ref = os.path.join(self.pele_dir, "ligand.pdb")
             self.native = cs.NATIVE.format(os.path.abspath(self.native), self.chain) if self.native else cs.NATIVE.format(os.path.abspath(self.ligand_ref), self.chain)
             self.glide_structs = os.path.join(self.pele_dir, "glide_calculations", "structures")
             self.topology = None if self.pdb else os.path.join(self.adap_ex_output, "topology.pdb")

        elif self.software in ["adaptive", "out_in", "induce_fit"]:
            self.adap_ex_output = None
            self.adap_ex_input = os.path.join(self.pele_dir, os.path.basename(self.system_fix))
            if self.software == "adaptive":
                self.ad_ex_temp = os.path.join(self.pele_dir, os.path.basename(self.adaptive))
                self.pele_exit_temp = os.path.join(self.pele_dir, os.path.basename(self.pele))
            else:
                self.ad_ex_temp = os.path.join(self.pele_dir, os.path.basename(self.files[0]))
                self.pele_exit_temp = os.path.join(self.pele_dir, os.path.basename(self.files[1]))
            self.ligand_ref = os.path.join(self.pele_dir, "ligand.pdb")
            self.native = cs.NATIVE.format(os.path.abspath(self.native), self.chain) if self.native else cs.NATIVE.format(os.path.abspath(self.ligand_ref), self.chain)
            self.topology = None if self.pdb else os.path.join("output", "topology.pdb")

    def create(self):
        if self.restart in cs.FIRST_RESTART:
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
        for file, destination_name in zip(self.files, self.file_names):
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
        if self.restart in ["all", "glide" ]:
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

