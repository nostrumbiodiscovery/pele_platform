import socket
import os

############################
# PRIVATE CONSTANTS
#############################

SCHRODINGER = os.environ.get("SCHRODINGER", "")
PELE = os.environ.get("PELE", "")

# DEFAULTS
# --------

COMPLEX = "complex.pdb"
RESULTS = "results"
LIG_RES = "LIG"
LIG_CHAIN = "Z"
FORCEFIELD = "OPLS2005"
PELE_CONFILE = "pele.conf"
CPUS = 140
RESTART = "true"
CLUSTERS = 40
PLATFORM_RESTART = "all"
EQ_STEPS = 50
GRIDRES = '10.0'
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
        
WATER_ENERGY =             '''
                            {{
                            "type": "bindingEnergy",\n\
                            "boundPartSelection": {{ "chains": {{ "names": ["{0}"] }} }},\n\
                            "tag": "water{0}"\n\
                            }},\n\
                           '''


UNBINDING = '''
            "modeMovingBox" : "unbinding",
            "exitCondition" : {{
                "type" : "metricMultipleTrajectories",
                "params" : {{
                       "metricCol" : {},
                       "exitValue" : {},
                       "condition" : "{}",
                       "numberTrajectories" : {}
            }}
            }},
            '''


DISTANCE_ATOMS =     '''
                     {{
                     "type":"com_distance",
                     "tag":"distance{2}",
                     "selection_group_1":{{
                     "atoms": {{ "ids":["{0}"]}}
                     }},
                     "selection_group_2":{{
                     "atoms": {{ "ids":["{1}"]}}
                     }}
                     }},
                     '''
BOX = '''

                "Box": {{
                    "type": "sphericalBox",
                    "radius": {0},
                    "fixedCenter": {1}
                }},
'''                     

WATER = '''
         "WaterPerturbation":
         {{
             "Box" :
             {{
                 "radius" : {},
                 "fixedCenter": {},
                 "type" : "sphericalBox"
             }},
             "watersToPerturb": {{ "links": {{ "ids": [ {} ] }} }},
             "parameters":
             {{
                 "temperature": {},
                 "numberOfStericTrials": {},
                 "overlapFactor": {},
                 "COMConstraintConstant": {}
             }}
         }}, 
'''


PCA = '''"preloadedModesIn" : "{}",'''


SELECTION_TO_PERTURB = '"selectionToPerturb" : { "chains" : { "names" : [ "$CHAIN" ] } },'
PERTURBATION = '''
          "Perturbation": {
                $BOX
                "perturbationType":"naive",
                "translationDirection": "steered",
                "rotationAngles": "nonCoupled",
                "parameters": {
                    "numberOfStericTrials": $STERIC_TRIALS,
                    "steeringUpdateFrequency": 0,
                    "overlapFactor": $OVERLAP
                }   
                
            },
'''
BE = '''
                        { "type": "bindingEnergy",

                           "boundPartSelection": { "chains": { "names": ["$CHAIN"] } }

                        },
'''

SASA='''
                        { "type": "sasa",

                           "tag": "sasaLig",

                           "selection": { "chains": { "names": ["$CHAIN"] } }

                        },
'''



LIGAND = '"ligandResname" : "$LIG_RES",'



#TEMPLATE KEYWORDS
#------------------

GLIDE_TEMPLATE = ["INPUT", "PRECISION"]

#RESTARTS:
#-----------

FIRST_RESTART = ["all",]
SECOND_RESTART = ["all", "adaptive"]
THIRD_RESTART = ["all", "adaptive", "pele"]
FOURTH_RESTART = ["all", "adaptive", "pele", "msm"] 

#PATHS
#-------

DIR = os.path.dirname(os.path.dirname(__file__))
ADAPTIVE_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "Adaptive/clusterAdaptiveRun.py"))

#MESSAGES&ERRORS
#-----------------

CLUSTER_ERROR = "Number of cpus ({}) must be bigger than clusters ({})"
SYSTEM = "\n\t**Missing residues found {}\n\t**Gaps found {}\n\t**Metals found {}"
