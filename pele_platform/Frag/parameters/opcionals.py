import pele_platform.Utilities.Helpers.constraints as cst
import pele_platform.constants.constants as cs

class FragOpcionalParameters():

    def __init__(self, args):
        #SIMULATION CONTROL
        self.frag_run = args.frag_run
        #CONSTRAINTS
        self.constraints = cst.retrieve_constraints(self.core, {}, {}, self.ca_constr)
        #CHAIN
        self.chain_core = args.chain_core if args.chain_core else self.simulation_params.get("chain_core", "L")
        #BOX
        self.box = cs.BOX.format(self.box_radius, self.box_center) if  self.box_radius else ""
        #LIGAND
        self.gridres = args.gridres
        self.plop_path = "PlopRotTemp_S_2017/ligand_prep.py"
        #OUTPUT
        self.criteria = args.frag_criteria if args.frag_criteria else self.simulation_params.get(
            "frag_criteria", "Binding Energy")
        self.output_folder = args.frag_output_folder if args.frag_output_folder else self.simulation_params.get(
            "frag_output_folder", "growing_steps")
        self.cluster_folder = args.frag_cluster_folder if args.frag_cluster_folder else self.simulation_params.get(
            "frag_cluster_folder", "clustering_PDBs")
        #OTHERS
        self.distcont = 4
        self.threshold = 0.3
        self.condition = "min"
        self.metricweights = "linear"
        self.nclusters = 5
        self.min_overlap = 0.5
        self.max_overlap = 0.7
        self.frag_chain = "L"
        self.banned = None
        self.limit = None
        self.rename = False
        self.threshold_clash = 1.7
        self.steering = 0
        self.translation_high = 0.05
        self.translation_low = 0.02
        self.rotation_high = 0.1
        self.rotation_low = 0.05
        self.explorative = False
        self.frag_radius = 10
        self.sampling_control = None
        self.only_prepare = False
        self.only_grow = False
        self.no_check = True
        

        
