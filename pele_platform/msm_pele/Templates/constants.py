import socket
import os

###########################
# USER TO CHANGE
###########################

machine = socket.getfqdn()
if "bsc.mn" in machine:
    SCHRODINGER = "/gpfs/projects/bsc72/SCHRODINGER_ACADEMIC"
    PELE = "/gpfs/projects/bsc72/MultiBoxPELE/bin/"
    PELE_BIN = "/gpfs/projects/bsc72/MultiBoxPELE/build/PELE-1.5"
    MPIRUN = "/apps/INTEL/2017.4/impi/2017.3.196/bin64"
    LICENSE = "/gpfs/projects/bsc72/PELE++/license"
    MMSHARE = None
    # Provisional workaround until best_struct.py is fixed
    ACCEPTED_STEPS_NAME = "numberOfAcceptedPeleSteps"
    CRITERIA = "sasaLig"

elif "mn.bsc" in machine:
    SCHRODINGER = "/gpfs/projects/bsc72/SCHRODINGER_ACADEMIC_NORD"
    PELE = "/gpfs/projects/bsc72/PELE++/nord/rev090518"
    PELE_BIN = "/gpfs/projects/bsc72/PELE++/nord/rev090518/bin/PELE-1.5_mpi"
    MPIRUN = "/apps/OPENMPI/1.8.1-mellanox/bin"
    LICENSE = "/gpfs/projects/bsc72/PELE++/license"
    MMSHARE = None
    # Provisional workaround until best_struct.py is fixed
    ACCEPTED_STEPS_NAME = "AcceptedSteps"
    CRITERIA = "SASA"

elif "bsccv" in machine:
    SCHRODINGER = "/data2/bsc72/SCHRODINGER_ACADEMIC"
    PELE = "/data/EAPM/PELE/PELE++/life/rev12489"
    PELE_BIN = "/data/EAPM/PELE/PELE++/life/rev12489/bin/PELE-1.5_mpi"
    MPIRUN = "/data2/apps/OPENMPI/1.6.1/bin"
    #LICENSE = "/gpfs/projects/bsc72/PELE++/license"
    MMSHARE = None
    LICENSE = "/data/EAPM/PELE/PELE++/license"
    # MMSHARE = "/data2/bsc72/SCHRODINGER_ACADEMIC/mmshare-v3.9/bin/Linux-x86_64"
    # SCHRODINGER_PYTHON_LIBS = "/data2/bsc72/SCHRODINGER_ACADEMIC/mmshare-v3.9/lib/Linux-x86_64/lib/python2.7/site-packages/"
    # Provisional workaround until best_struct.py is fixed
    ACCEPTED_STEPS_NAME = "AcceptedSteps"
    CRITERIA = "SASA"
elif "NBD" in machine:
    SCHRODINGER = "/sNow/easybuild/centos/7.4.1708/Skylake/software/schrodinger-2017/"
    PELE = "/sNow/easybuild/centos/7.4.1708/Skylake/software/PELE/1.5.0.2524/"
    PELE_BIN = "/home/dsoler/cleanPELE_rev/build_gnu/Pele_mpi"
    MPIRUN = "/sNow/easybuild/centos/7.4.1708/Skylake/software/OpenMPI/2.1.2-GCC-6.4.0-2.28/bin/"
    LICENSE = "/sNow/easybuild/centos/7.4.1708/Skylake/software/PELE/licenses/"
    MMSHARE = None
    # Provisional workaround until best_struct.py is fixed
    ACCEPTED_STEPS_NAME = "numberOfAcceptedPeleSteps"
    CRITERIA = "sasaLig"
else:
    SCHRODINGER = "$SCHRODINGER"
    PELE = "$PELE"
    PELE_BIN = "$PELE_BIN"
    MPIRUN = "$MPIRUN"
    LICENSE = "$LICENSE"
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
ADAPTIVE_KEYWORDS = ["RESTART", "OUTPUT", "INPUT", "CPUS", "PELE_CFILE", "LIG_RES", "SEED", "STEPS", "ITERATIONS", "MSM_CLUST", "LAGTIME", "MIN_POS"]
EX_ADAPTIVE_KEYWORDS = ["RESTART", "OUTPUT", "INPUT", "CPUS", "PELE_CFILE", "LIG_RES", "EQ_STEPS", "SEED", "EXIT_ITERS", "EQ_STRUCT"]
EX_PELE_KEYWORDS = ["NATIVE", "FORCEFIELD", "CHAIN", "CONSTRAINTS", "LICENSES", "LOGFILE", "SOLVENT"]
PELE_KEYWORDS = [ "RESTART", "OUTPUT", "INPUT", "SEED", "STEPS", "BOX", "BOX_METRIC", "SASA_min", "SASA_max", "TEMP" ]
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
BOX_METRIC = '''
            {
                "type": "isPerturbedAtomSetCOMOutOfTheBox"
            }
         '''

BOX = '''
MODEL $MODEL
HEADER    CORNERS OF BOX
REMARK    CENTER (X Y Z)   $CENTER_X  $CENTER_Y  $CENTER_Z
REMARK    RADIUS $RADIUS
REMARK    DIMENSIONS (X Y Z)   31.116  26.557  29.958
ATOM      1  DUA BOX     1 $V1
ATOM      2  DUB BOX     1 $V2
ATOM      3  DUC BOX     1 $V3
ATOM      4  DUD BOX     1 $V4
ATOM      5  DUE BOX     1 $V5
ATOM      6  DUF BOX     1 $V6
ATOM      7  DUG BOX     1 $V7
ATOM      8  DUH BOX     1 $V8
CONECT    1    2    4    5
CONECT    2    1    3    6
CONECT    3    2    4    7
CONECT    4    1    3    8
CONECT    5    1    6    8
CONECT    6    2    5    7
CONECT    7    3    6    8
CONECT    8    4    5    7
ENDMDL
'''

INPUT_PELE = '{{ "files" : [ {{ "path" : "{}" }} ] }}'


SYSTEM = "System {} checked successfully\n\t**Missing residues found {}\n\t**Gaps found {}\n\t**Metals found {}"



# FOLDERS&PATH
DIR = os.path.dirname(__file__)
ADAPTIVE_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "Adaptive/clusterAdaptiveRun.py"))
FOLDERS = ["",
           "DataLocal/Templates/OPLS2005/HeteroAtoms/",
           "DataLocal/Templates/AMBER99sb/HeteroAtoms/",
           "DataLocal/Templates/AMBER99sbBSC0/HeteroAtoms/",
           "DataLocal/LigandRotamerLibs",
           "DataLocal/OBC",
           "output_pele",
           "output_adaptive_exit",
           "output_clustering",
           "results"
          ]

FILES_SP = [os.path.join(DIR, "Templates/box.pdb"), os.path.join(DIR, "Templates/pele_SP.conf"),
                 os.path.join(DIR, "Templates/adaptive_exit.conf"),
                 os.path.join(DIR, "Templates/pele_exit.conf")]

FILES_XP = [os.path.join(DIR, "Templates/box.pdb"), os.path.join(DIR, "Templates/pele_XP.conf"),
                 os.path.join(DIR, "Templates/adaptive_exit.conf"),
                 os.path.join(DIR, "Templates/pele_exit.conf")]

FILES_XP2 = [os.path.join(DIR, "Templates/box.pdb"), os.path.join(DIR, "Templates/pele_XP2.conf"),
                 os.path.join(DIR, "Templates/adaptive_exit.conf"),
                os.path.join(DIR, "Templates/pele_exit.conf")]

FILES_TEST_XP = [os.path.join(DIR, "Templates/box.pdb"), os.path.join(DIR, "Templates/pele_XP.conf"),
                 os.path.join(DIR, "Templates/adaptive_exit_test.conf"),
                 os.path.join(DIR, "Templates/pele_exit.conf")]

FILES_TEST = [os.path.join(DIR, "Templates/box.pdb"), os.path.join(DIR, "Templates/pele_SP.conf"),
                 os.path.join(DIR, "Templates/adaptive_exit_test.conf"),
                 os.path.join(DIR, "Templates/pele_exit.conf")]

FILES_NAME = ["box.pdb", "pele.conf", "adaptive_exit.conf",  "pele_exit.conf"]

# ERRORS
CLUSTER_ERROR = "Number of cpus ({}) must be bigger than clusters ({})"

