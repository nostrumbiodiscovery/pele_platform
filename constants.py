import socket
import os

###########################
# USER TO CHANGE
###########################

machine = socket.getfqdn()
if "bsc.mn" in machine:
    SCHRODINGER = "/gpfs/projects/bsc72/SCHRODINGER_ACADEMIC"
    PELE = "/gpfs/projects/bsc72/PELE++/mniv/rev12455"
    ADAPTIVE = "/gpfs/projects/bsc72/adaptiveSampling/bin/v1.4.2"
    MPIRUN = "/apps/INTEL/2017.4/impi/2017.3.196/bin64"
    LICENSE = "/gpfs/projects/bsc72/PELE++/license"
    MMSHARE = None
    # Provisional workaround until best_struct.py is fixed
    ACCEPTED_STEPS_NAME = "numberOfAcceptedPeleSteps"
    CRITERIA = "sasaLig"

elif "mn.bsc" in machine:
    SCHRODINGER = "/gpfs/projects/bsc72/SCHRODINGER_ACADEMIC_NORD"
    PELE = "/gpfs/projects/bsc72/PELE++/nord/rev12489"
    ADAPTIVE = "/gpfs/projects/bsc72/adaptiveSampling/bin_nord/v1.4.2_schr"
    MPIRUN = "/apps/OPENMPI/1.8.1-mellanox/bin"
    LICENSE = "/gpfs/projects/bsc72/PELE++/license"
    MMSHARE = None
    # Provisional workaround until best_struct.py is fixed
    ACCEPTED_STEPS_NAME = "AcceptedSteps"
    CRITERIA = "SASA"

elif "bsccv" in machine:
    SCHRODINGER = "/data2/bsc72/SCHRODINGER_ACADEMIC"
    PELE = "/data/EAPM/PELE/PELE++/life/rev12489"
    ADAPTIVE = "/data2/bsc72/AdaptiveSampling/bin/v1.4.2_bis"
    MPIRUN = "/data2/apps/OPENMPI/1.6.1/bin"
    #LICENSE = "/gpfs/projects/bsc72/PELE++/license"
    MMSHARE = None
    LICENSE = "/data/EAPM/PELE/PELE++/license"
    # MMSHARE = "/data2/bsc72/SCHRODINGER_ACADEMIC/mmshare-v3.9/bin/Linux-x86_64"
    # SCHRODINGER_PYTHON_LIBS = "/data2/bsc72/SCHRODINGER_ACADEMIC/mmshare-v3.9/lib/Linux-x86_64/lib/python2.7/site-packages/"
    # Provisional workaround until best_struct.py is fixed
    ACCEPTED_STEPS_NAME = "AcceptedSteps"
    CRITERIA = "SASA"
else:
    SCHRODINGER = "/sNow/easybuild/centos/7.4.1708/Skylake/software/schrodinger2017-4/"
    PELE = "/sNow/easybuild/centos/7.4.1708/Skylake/software/PELE/1.5.0.2524/"
    ADAPTIVE = "/home/dsoler/repos/AdaptivePELE/"
    MPIRUN = "/sNow/easybuild/centos/7.4.1708/Skylake/software/OpenMPI/1.8.4-GCC-4.9.2/bin"
    LICENSE = "/sNow/easybuild/centos/7.4.1708/Skylake/software/PELE/licenses/"
    MMSHARE = None
    # Provisional workaround until best_struct.py is fixed
    ACCEPTED_STEPS_NAME = "numberOfAcceptedPeleSteps"
    CRITERIA = "sasaLig"



############################
# PRIVATE CONSTANTS
#############################

# DEFAULT VALUES
COMPLEX = "complex.pdb"
RESULTS = "results"
LIG_RES = "LIG"
LIG_CHAIN = "Z"
FORCEFIELD = "OPLS2005"
PELE_CONFILE = "pele.conf"
CPUS = 140
RESTART = True
CLUSTERS = 40
PLATFORM_RESTART = "all"
EQ_STEPS = 50
GRIDRES = '10.0'

# TEMPLATE KEYWORDS
ADAPTIVE_KEYWORDS = ["RESTART", "OUTPUT", "INPUT", "CPUS", "PELE_CFILE", "LIG_RES", "SEED"]
EX_ADAPTIVE_KEYWORDS = ["RESTART", "OUTPUT", "INPUT", "CPUS", "PELE_CFILE", "LIG_RES", "EQ_STEPS", "SEED"]
EX_PELE_KEYWORDS = ["NATIVE", "FORCEFIELD", "CHAIN", "CONSTRAINTS", "CPUS", "LICENSES"]
PELE_KEYWORDS = ["BOX_CENTER", "BOX_RADIUS"]
NATIVE = '''
                        {{

                           "type": "rmsd",

                           "Native": {{\n\
                            "path":\n\
                            "{}" }},\n\
    
                           "selection": {{ "chains": {{ "names": [ "{}" ] }} }},\n\

                           "includeHydrogens": false,\n\

                           "doSuperposition": false,\n\

                           "tag" : "ligandRMSD"\n\

                        }},\n\


'''

SYSTEM = "System {} checked successfully\n\t**Missing residues found {}\n\t**Gaps found {}\n\t**Metals found {}"



# FOLDERS&PATH
DIR = os.path.dirname(__file__)
ADAPTIVE_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "Adaptive/clusterAdaptiveRun.py"))
FOLDERS = ["",
           "DataLocal/Templates/OPLS2005/HeteroAtoms/",
           "DataLocal/Templates/AMBER99sb/HeteroAtoms/",
           "DataLocal/Templates/AMBER99sbBSC0/HeteroAtoms/",
           "DataLocal/LigandRotamerLibs",
           "output_pele",
           "output_adaptive_exit",
           "output_clustering"
          ]

FILES_SP = [os.path.join(DIR, "Templates/box.pdb"), os.path.join(DIR, "Templates/pele_SP.conf"),
                 os.path.join(DIR, "Templates/adaptive_exit.conf"), os.path.join(DIR, "Templates/adaptive_long.conf"),
                 os.path.join(DIR, "Templates/pele_exit.conf")]

FILES_XP = [os.path.join(DIR, "Templates/box.pdb"), os.path.join(DIR, "Templates/pele_XP.conf"),
                 os.path.join(DIR, "Templates/adaptive_exit.conf"), os.path.join(DIR, "Templates/adaptive_long.conf"),                                                                                             
                 os.path.join(DIR, "Templates/pele_exit.conf")]

FILES_TEST = [os.path.join(DIR, "Templates/box.pdb"), os.path.join(DIR, "Templates/pele_SP.conf"),
                 os.path.join(DIR, "Templates/adaptive_exit_test.conf"), os.path.join(DIR, "Templates/adaptive_long_test.conf"),
                 os.path.join(DIR, "Templates/pele_exit.conf")]

FILES_NAME = ["box.pdb", "pele.conf", "adaptive_exit.conf", "adaptive_long.conf", "pele_exit.conf"] 

# ERRORS
CLUSTER_ERROR = "Number of cpus ({}) must be bigger than clusters ({})"
