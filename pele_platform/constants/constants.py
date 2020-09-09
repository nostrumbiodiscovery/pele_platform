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
             "watersToPerturb": {{ "links": {{ "ids": [ {} ] }} }},
             "parameters":
             {{
                 {}
                 "temperature": {},
                 "numberOfStericTrials": {},
                 "overlapFactor": {},
                 "COMConstraintConstant": {}
             }},
             "waterSites": {}
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
                    "steeringUpdateFrequency": $STEERING,
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

# PPI waters
water_O =  "HETATM {}  OW  HOH {}{:>4}     {}  1.00  0.00           O\n"
water_H1 = "HETATM {}  1HW HOH {}{:>4}     {}  1.00  0.00           H\n"
water_H2 = "HETATM {}  2HW HOH {}{:>4}     {}  1.00  0.00           H\n"
water = [water_O, water_H1, water_H2]


# METALS
metals = ['LI', 'BE', 'NA', 'MG', 'AL', 'K', 'CA', 'SC', 'TI', 'V', 'CR', 'MN', 'FE', 'CO', 'NI', 'CU', 'ZN', 'GA', 'RB', 'SR', 'Y', 'ZR', 'NB', 'MO', 'TC', 'RU', 'RH', 'PD', 'AG', 'CD', 'IN', 'SN', 'CS', 'BA', 'LA', 'CE', 'PR', 'ND', 'PM', 'SM', 'EU', 'GD', 'TB', 'DY', 'HO', 'ER', 'TM', 'YB', 'LU', 'HF', 'TA', 'W', 'RE', 'OS', 'IR', 'PT', 'AU', 'HG', 'TL', 'PB', 'BI', 'PO', 'FR', 'RA', 'AC', 'TH', 'PA', 'U', 'NP', 'PU', 'AM', 'CM', 'BK', 'CF', 'ES', 'FM', 'MD', 'NO', 'LR', 'RF', 'DB', 'SG', 'BH', 'HS', 'MT', 'DS', 'RG', 'CN', 'NH', 'FL', 'MC', 'LV']


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
