"""
This file contains de perturbation
parameters of PELE for each type of simulation
"""



BIAS = ''' 
             ,
             "parametersChanges" : [
                 { "ifAnyIsTrue": [ "rand >= .5" ],
                     "doThesechanges": { "Perturbation::parameters": { "rotationScalingFactor": 0.1 } },
                     "otherwise": { "Perturbation::parameters": { "rotationScalingFactor": 0.25 } }
                 },
                 { "ifAnyIsTrue": [ "rand1 >= 0.5" ],
                     "doThesechanges": { "Perturbation::parameters": { "translationRange": 3.0, "numberOfTrials" : 10 } },
                     "otherwise": { "Perturbation::parameters": { "translationRange": 1.25, "numberOfTrials" : 10 } }
                 },
                 { "ifAnyIsTrue": [ "sasaLig >= 0.85" ],
                     "doThesechanges": { "Perturbation::parameters": { "translationRange": 5.0, "numberOfTrials" : 20 } },
                     "otherwise": {  }
                 },
                 { "ifAnyIsTrue": [ "sasaLig <= 0.35" ],
                     "doThesechanges": { "Perturbation::parameters": { "rotationScalingFactor": 0.1, "translationRange": 0.75, "numberOfTrials" : 40, "steeringUpdateFrequency": 0} },
                     "otherwise": {  }
                 }
]
'''

IN_OUT = '''
             ,

                 "parametersChanges" : [

                     { "ifAnyIsTrue": [ "rand >= .5" ],

                         "doThesechanges": { "Perturbation::parameters": { "rotationScalingFactor": 0.05 } },

                         "otherwise": { "Perturbation::parameters": { "rotationScalingFactor": 0.15 } }

                     },

                     {

                       "ifAnyIsTrue": [ "rand1 >= 0.40" ],

                         "doThesechanges": { "Perturbation::parameters": { "translationRange": 2.0 } },

                         "otherwise": { "Perturbation::parameters": { "translationRange": 2.0 } }

                    },
                    {
                      "ifAnyIsTrue": [ "rand3 >= 0.10" ],
                        "doThesechanges": { "Perturbation::parameters": { "steeringUpdateFrequency": 1, "numberOfTrials" : 10 } },
                        "otherwise": { "Perturbation::parameters": { "steeringUpdateFrequency": 0, "numberOfTrials" : 25 } }
                    },
                    {
                      "ifAnyIsTrue": [ "sasaLig <= 0.15" ],
                         "doThesechanges": { "Pele::parameters": { "perturbationCOMConstraintConstant" : 0.25 }, "Perturbation::parameters": { "translationRange": 1.0 }},
                         "otherwise": { "Pele::parameters": { "perturbationCOMConstraintConstant" : 1.0 } }
                    }, 
                    {
                      "ifAnyIsTrue": [ "sasaLig >= 0.75" ],
                         "doThesechanges": { "Pele::parameters": { "perturbationCOMConstraintConstant" : 10.0 }, "Perturbation::parameters": { "steeringUpdateFrequency": 1, "numberOfTrials" : 4 }},
                         "otherwise": { }
                    } 

                  ]

'''

RESCORING = '''
             ,
             "parametersChanges" : [
             
                  { "ifAnyIsTrue": [ "rand >= .5" ],
             
                      "doThesechanges": { "Perturbation::parameters": { "rotationScalingFactor": 0.1 } },
             
                      "otherwise": { "Perturbation::parameters": { "rotationScalingFactor": 0.05 } }
             
                  },
             
                  { "ifAnyIsTrue": [ "rand1 >= 0.5" ],
             
                      "doThesechanges": { "Perturbation::parameters": { "translationRange": 0.25} },
             
                      "otherwise": { "Perturbation::parameters": { "translationRange": 0.5} }
             
                  },
             
                  {  "ifAnyIsTrue": [ "rand2 >= 0.5" ],
             
                         "doThesechanges": {  "Perturbation::parameters": { "steeringUpdateFrequency": 0, "numberOfTrials": 25 } },
             
                         "otherwise": {  "Perturbation::parameters": { "steeringUpdateFrequency": 0 , "numberOfTrials": 25  }}
             
                 }
             
             ]
'''



INDUCED_FIT = '''
             ,
"parametersChanges" : [
     { "ifAnyIsTrue": [ "rand >= .5" ],
         "doThesechanges": { "Perturbation::parameters": { "rotationScalingFactor": 0.1 } },
         "otherwise": { "Perturbation::parameters": { "rotationScalingFactor": 0.25 } }
     },
     { "ifAnyIsTrue": [ "rand1 >= 0.5" ],
         "doThesechanges": { "Perturbation::parameters": { "translationRange": 0.5} },
         "otherwise": { "Perturbation::parameters": { "translationRange": 1.0} }
     },
     {  "ifAnyIsTrue": [ "rand2 >= 0.5" ],
            "doThesechanges": {  "Perturbation::parameters": { "steeringUpdateFrequency": 0, "numberOfTrials": 30 } },
            "otherwise": {  "Perturbation::parameters": { "steeringUpdateFrequency": 1 , "numberOfTrials": 10  }}
    }
]
'''


GLOBAL = '''
             ,
             "parametersChanges" : [
                 { "ifAnyIsTrue": [ "rand >= .5" ],
                     "doThesechanges": { "Perturbation::parameters": { "rotationScalingFactor": 0.1 } },
                     "otherwise": { "Perturbation::parameters": { "rotationScalingFactor": 0.25 } }
                 },
                 { "ifAnyIsTrue": [ "rand1 >= 0.5" ],
                     "doThesechanges": { "Perturbation::parameters": { "translationRange": 3.0, "numberOfTrials" : 10 } },
                     "otherwise": { "Perturbation::parameters": { "translationRange": 1.25, "numberOfTrials" : 10 } }
                 },
                 { "ifAnyIsTrue": [ "sasaLig >= 0.85" ],
                     "doThesechanges": { "Perturbation::parameters": { "translationRange": 5.0, "numberOfTrials" : 20 } },
                     "otherwise": {  }
                 },
                 { "ifAnyIsTrue": [ "sasaLig <= 0.35" ],
                     "doThesechanges": { "Perturbation::parameters": { "rotationScalingFactor": 0.1, "translationRange": 0.75, "numberOfTrials" : 40, "steeringUpdateFrequency": 0} },
                     "otherwise": {  }
                 }
              ]
'''

WATER_BS = ""

WATER_LIG = '''
                 ,
                 "parametersChanges" :
                 [
                     { "ifAnyIsTrue": [ "rand1 >= 0.5" ],
                         "doThesechanges": { "Perturbation::parameters": { "rotationScalingFactor": 0.01 } },
                         "otherwise": { "Perturbation::parameters": { "rotationScalingFactor": 0.025 } }
                     },
                     {
                       "ifAnyIsTrue": [ "rand2 >= 0.5" ],
                         "doThesechanges": { "Perturbation::parameters": { "translationRange": 0.75 } },
                         "otherwise": { "Perturbation::parameters": { "translationRange": 0.5} }
                     },
                     { "ifAnyIsTrue": [ "rand3 <= 0.5 " ],
                         "doThesechanges": { "WaterPerturbation::parameters": { "translationRange": 3, "temperature": 5000 } },
                         "otherwise": { "WaterPerturbation::parameters": { "translationRange": 4, "temperature": 7500 } }
                     },
                     { "ifAnyIsTrue": [ "rand3 > 0.85 " ],
                         "doThesechanges": { "WaterPerturbation::parameters": { "translationRange": 5, "temperature": 10000 } },
                         "otherwise": {}
                     }
                  ]
'''

WATER_LIG_EXPL = '''
                 ,
                 "parametersChanges" : [
                     { "ifAnyIsTrue": [ "rand1 >= 0.5" ],
                         "doThesechanges": { "Perturbation::parameters": { "rotationScalingFactor": 0.1 } },
                         "otherwise": { "Perturbation::parameters": { "rotationScalingFactor": 0.2 } }
                     },
                     {
                       "ifAnyIsTrue": [ "rand2 >= 0.5" ],
                         "doThesechanges": { "Perturbation::parameters": { "translationRange": 1.75 } },
                         "otherwise": { "Perturbation::parameters": { "translationRange": 0.75} }
                     },
                     {
                         "ifAnyIsTrue": [ "rand3 <= 0.4" ],
                         "doThesechanges": { "WaterPerturbation::parameters": { "translationRange": 4.0 } },
                         "otherwise": { "WaterPerturbation::parameters": { "translationRange": 3.0} }
                     },
                    {
                      "ifAnyIsTrue": [ "rand4 >= 0.50" ],
                        "doThesechanges": { "Perturbation::parameters": { "steeringUpdateFrequency": 1, "numberOfTrials" : 10 } },
                        "otherwise": { "Perturbation::parameters": { "steeringUpdateFrequency": 0, "numberOfTrials" : 25 } }
                    },
                     {
                         "ifAnyIsTrue": [ "rand3 >= 0.85" ],
                         "doThesechanges": { "WaterPerturbation::parameters": { "translationRange": 2.0 } }
                     }
                  ]
'''

GLIDE = '''
            ,
			"parametersChanges" : [
				{ "ifAnyIsTrue": [ "rand >= 0.5" ], 
					"doThesechanges": { 
						"Perturbation::parameters": { "rotationScalingFactor": 0.30 } 
					}, 
					"otherwise": { 
						"Perturbation::parameters": { "rotationScalingFactor": 0.05 } 
					} 
				}, 
				{ "ifAnyIsTrue": [ "rand1 <= 0.20" ], 
					"doThesechanges": { 
						"Perturbation::parameters": { "translationRange": 0.5 } 
					}, 
					"otherwise": { "Perturbation::parameters": { "translationRange": 1.5 } } 
				}, 
                                { "ifAnyIsTrue": [ "rand2 <= 0.5 " ], 
					"doThesechanges": { 
						"Perturbation::parameters": { 
							"steeringUpdateFrequency": 0,
							"numberOfTrials": 10}
					}, 
					"otherwise": { 
						"Perturbation::parameters": { 
							"steeringUpdateFrequency": 1, 
							"numberOfTrials": 10}
					} 
				}, 
				{ "ifAnyIsTrue": [ "hbond_hinge <= 3.0 " ], 
					"doThesechanges": { 
						"Perturbation::parameters": { 
							"translationRange": 0.5, 
							"steeringUpdateFrequency": 0, 
							"numberOfTrials": 25 }
					},
                    "otherwise": {} 
				} 
			] 
'''

FRAG = '''
,
                 "parametersChanges" : [
                     { "ifAnyIsTrue": [ "rand >= 0.5" ],
                         "doThesechanges": { "Perturbation::parameters": { "rotationScalingFactor": $ROTATION_HIGH } },
                         "otherwise": { "Perturbation::parameters": { "rotationScalingFactor": $ROTATION_LOW } }
                     },
                     { "ifAnyIsTrue": [ "rand1 <= 0.5 " ],
                         "doThesechanges": { "Perturbation::parameters": { "translationRange": $TRANSLATION_LOW } },
                         "otherwise": { "Perturbation::parameters": { "translationRange": $TRANSLATION_HIGH } }
                     }
                  ]
'''
