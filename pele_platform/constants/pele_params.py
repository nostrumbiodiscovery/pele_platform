



BIAS = ''' 
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
