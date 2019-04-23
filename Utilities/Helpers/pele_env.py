import os
import shutil
import warnings
import random
import logging
import PELE_Platform.constants as cs
import PELE_Platform.features as fs
import PELE_Platform.Utilities.Helpers.helpers as hp



class EnviroBuilder(object):
    """
        Base class wher the needed pele environment
        is build by creating folders and files
    """

    def __init__(self, args):
        self.build_variables(args)
        self.build_paths(args)

    @classmethod
    def build_env(cls, args):
        env = cls(args)
        env.create()
        return env

    def build_variables(self, args):

        self.system = args.system
        self.forcefield = args.forcefield
        self.residue = args.residue
        self.skip_prep = args.skip_prep
        self.templates = os.path.abspath(os.path.join(os.path.dirname(os.path.dirname(__file__)), "PeleTemplates"))
        self.solvent = args.solvent
        self.cpus = args.cpus = args.cpus if not args.test else 4
        self.adapt_conf = args.adapt_conf
        self.confile = args.confile
        self.input = args.input
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
        self.nonstandard = args.nonstandard
        self.external_template = args.template
        self.external_rotamers = args.rotamers
        self.lagtime = args.lagtime
        self.msm_clust = args.msm_clust
        self.seed = random.randrange(1, 70000)
        self.license = '''"{}"'''.format(cs.LICENSE)
        self.equil_steps = 1 if self.test else int(cs.EQ_STEPS/self.cpus) + 1 #+1 to avoid being 0
        self.sasa_max = None
        self.sasa_min = None
        self.equilibration = "true" if args.noeq else "false"
        self.box_radius = None
        self.box_center = "["+ ",".join(args.user_center) + "]" if args.user_center else None
        self.constraints = None
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
            

    def build_paths(self, args):


        pele_dir = os.path.abspath("{}_Pele".format(self.residue))

        if not self.folder:
            self.pele_dir = hp.is_repited(pele_dir) if self.restart in cs.FIRST_RESTART else hp.is_last(pele_dir)
        else:
            self.pele_dir = os.path.abspath(self.folder)

        if self.mae_lig:
            self.system_fix = os.path.join(self.pele_dir, "{}_complex_processed.pdb".format(os.path.splitext(os.path.basename(self.system))[0]))
        else:
            self.system_fix = os.path.join(self.pele_dir, "{}_processed.pdb".format(os.path.splitext(os.path.basename(self.system))[0]))

        self.rotamers_folder = os.path.join(self.pele_dir, "DataLocal/LigandRotamerLibs/")
        self.template_folder = os.path.join(self.pele_dir, "DataLocal/Templates/{}/HeteroAtoms/".format(self.forcefield))


        self.adap_ex_input = os.path.join(self.pele_dir, os.path.basename(self.system_fix))
        self.obc_tmp = os.path.join(cs.DIR, "Templates/solventParamsHCTOBC.txt")
        self.obc_file = os.path.join(self.pele_dir, "DataLocal/OBC/solventParamsHCTOBC.txt")
        self.cluster_output = os.path.join(self.pele_dir, "output_clustering")
        self.clusters_output = os.path.join(self.cluster_output, "clusters_{}_KMeans_allSnapshots.pdb".format(self.clusters))
        self.receptor = os.path.join(self.pele_dir, "receptor.pdb")
        self.ligand_ref = os.path.join(self.pele_dir, "ligand.pdb")
        self.topology = None if self.pdb else os.path.join("output_pele", "topology.pdb")
        self.native = cs.NATIVE.format(os.path.abspath(self.native), self.chain) if self.native else ""
        self.pele_temp = os.path.join(self.pele_dir, "pele.conf")
        self.adap_l_input = "{}/initial_*"
        self.adap_l_output = os.path.join(self.pele_dir, "output_pele")
        self.ad_l_temp = os.path.join(self.pele_dir, "adaptive_long.conf")
        self.box_temp = os.path.join(self.pele_dir, "box.pdb")
        self.glide_template = os.path.join(self.pele_dir, "glide.in") 
        self.glide_structs = os.path.join(self.pele_dir, "glide_calculations", "structures")

        #####Define files and folders HIDING VARIABLES TO CHANGE#########
        for key, value in fs.retrieve_software_settings(args, self.pele_dir).items():
            setattr(self, key, value)


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
