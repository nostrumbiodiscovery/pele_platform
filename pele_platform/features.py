import os
import pele_platform.constants.constants as cs
import pele_platform.constants.pele_params as pcs

"""
Description of the file: This file specifies what files to use
for each job.

"""
    
def retrieve_software_settings(args, pele_dir):

        SOFTWARE_CONSTANTS = {
             "msm": {
                 "adap_ex_output" : os.path.join(pele_dir, "output_adaptive_exit",),
                 "ad_ex_temp" : os.path.join(pele_dir, "adaptive_exit.conf"),
                 "pele_exit_temp" : os.path.join(pele_dir, "pele_exit.conf"),
                 "folders" : ["",
                           "DataLocal/Templates/OPLS2005/HeteroAtoms/",
                           "DataLocal/Templates/AMBER99sb/HeteroAtoms/",
                           "DataLocal/Templates/AMBER99sbBSC0/HeteroAtoms/",
                           "DataLocal/LigandRotamerLibs",
                           "DataLocal/OBC",
                           "output_pele",
                           "output_adaptive_exit",
                           "output_clustering"
                          ], 
                 "file_names" :  ["box.pdb", "pele.conf", "adaptive_exit.conf", 
                                  "adaptive_long.conf", "pele_exit.conf"],

                 "files" : { "SP" :  [os.path.join(cs.DIR, "Templates/box.pdb"), os.path.join(cs.DIR, "Templates/pele_SP.conf"),
                                 os.path.join(cs.DIR, "Templates/adaptive_exit.conf"),
                                 os.path.join(cs.DIR, "Templates/adaptive_long.conf"),
                                 os.path.join(cs.DIR, "Templates/pele_exit.conf")],
                 
                             "XP" :   [os.path.join(cs.DIR, "Templates/box.pdb"), os.path.join(cs.DIR, "Templates/pele_XP.conf"),
                                 os.path.join(cs.DIR, "Templates/adaptive_exit.conf"),
                                 os.path.join(cs.DIR, "Templates/adaptive_long.conf"),
                                 os.path.join(cs.DIR, "Templates/pele_exit.conf")],
                             "test" : [os.path.join(cs.DIR, "Templates/box.pdb"), os.path.join(cs.DIR, "Templates/pele_SP.conf"),
                                 os.path.join(cs.DIR, "Templates/adaptive_exit_test.conf"),
                                 os.path.join(cs.DIR, "Templates/adaptive_long_test.conf"),
                                 os.path.join(cs.DIR, "Templates/pele_exit.conf")]
                                 },
                 "simulation_params" : {}
                 },
            
             "glide": {
                 "adap_ex_output" : os.path.join( pele_dir, "output_adaptive"),
                 "ad_ex_temp" : os.path.join( pele_dir, "adaptive.conf"),
                 "pele_exit_temp" : os.path.join( pele_dir, "pele.conf"),
                 "folders": [ "",
                    "DataLocal/Templates/OPLS2005/HeteroAtoms/",
                    "DataLocal/Templates/AMBER99sb/HeteroAtoms/",
                    "DataLocal/Templates/AMBER99sbBSC0/HeteroAtoms/",
                    "DataLocal/LigandRotamerLibs",
                    "DataLocal/OBC",
                    "output_adaptive",
                    "output_clustering",
                    "glide_calculations/structures"
                    ],
                  "file_names" : ["glide.in", "adaptive.conf", "pele.conf"],
                  "files" : { "glide" : [os.path.join(cs.DIR, "Templates/glide.in"),
                                         os.path.join(cs.DIR, "Templates/template_adaptive.conf"),
                                         os.path.join(cs.DIR, "Templates/pele_template.conf") ]
                            },
                 "simulation_params" : {
                             "glide": {"spawning_type": "epsilon", "bias_column": 6, "epsilon":0.75, "density": "null",
                                      "simulation_type": "pele", "iterations": 30, "pele_steps": 12, 
                                      "cluster_values": "[2, 3, 4, 5]", "cluster_conditions": "[1.25, 0.75, 0.5]",
                                      "steric_trials": 500, "overlap_factor": 0.65, "params": pcs.GLIDE,
                                      "box_radius": 10},
                            }

                 },
            
             "adaptive" :  {
                 "adap_ex_output" : None,
                 "ad_ex_temp" : os.path.join( pele_dir, "adaptive.conf"),
                 "pele_exit_temp" : os.path.join( pele_dir, "pele.conf"),
                 "folders" : [ "",
                    "DataLocal/Templates/OPLS2005/HeteroAtoms/",
                    "DataLocal/Templates/AMBER99sb/HeteroAtoms/",
                    "DataLocal/Templates/AMBER99sbBSC0/HeteroAtoms/",
                    "DataLocal/LigandRotamerLibs",
                    "DataLocal/OBC"
                    ],
                  "file_names" : ["adaptive.conf", "pele.conf" ],
                  "files" : { "adaptive" : [ args.adaptive, args.pele ],
                              "full" : [ os.path.join(cs.DIR, "Templates/template_adaptive.conf"), 
                                           os.path.join(cs.DIR, "Templates/pele_template.conf") ],
                              "induce_fit" : [ os.path.join(cs.DIR, "Templates/template_adaptive.conf"),
                                           os.path.join(cs.DIR, "Templates/pele_template.conf") ],
                              "in_out" : [ os.path.join(cs.DIR, "Templates/template_adaptive.conf"),
                                           os.path.join(cs.DIR, "Templates/pele_template.conf") ],
                              "in_out_soft" : [ os.path.join(cs.DIR, "Templates/template_adaptive.conf"),
                                           os.path.join(cs.DIR, "Templates/pele_template.conf") ],
                              "water_exp":  [ os.path.join(cs.DIR, "Templates/template_adaptive.conf"),
                                           os.path.join(cs.DIR, "Templates/pele_template.conf") ],
                              "water_lig":  [ os.path.join(cs.DIR, "Templates/template_adaptive.conf"),
                                           os.path.join(cs.DIR, "Templates/pele_template.conf") ],
                              "bias":  [ os.path.join(cs.DIR, "Templates/template_adaptive.conf"),
                                           os.path.join(cs.DIR, "Templates/pele_template.conf") ]
                            },
                 "simulation_params" : {
                             "adaptive": {"spawning_type": "epsilon", "bias_column": 5, "epsilon":0.15, "density": "continuous",
                                      "simulation_type": "pele", "iterations": 100, "pele_steps": 4, 
                                      "cluster_values": "[2.5, 5, 7]", "cluster_conditions": "[1, 0.6, 0.0]",
                                      "steric_trials": 500, "overlap_factor": 0.65, "params": pcs.BIAS,
                                      "box_radius": 10},
                             "full": {"spawning_type": "inverselyProportional", "bias_column": 5, "epsilon":0.25, "density": "continuous",
                                      "simulation_type": "pele", "iterations": 100, "pele_steps": 4, 
                                      "cluster_values": "[2.5, 5, 7]", "cluster_conditions": "[1, 0.6, 0.0]",
                                      "steric_trials": 200, "overlap_factor": 0.65, "params": pcs.GLOBAL,
                                      "box_radius":None},
                             "induce_fit" : {"spawning_type": "inverselyProportional", "bias_column": 5, "epsilon":0.25, "density": "null",
                                      "simulation_type": "pele", "iterations": 50, "pele_steps": 12, 
                                      "cluster_values": "[1.75, 2.5, 4, 6]", "cluster_conditions": "[1, 0.6, 0.4, 0.0]",
                                      "steric_trials": 500, "overlap_factor": 0.65, "params": pcs.INDUCED_FIT,
                                      "box_radius": 6},
                             "in_out" : {"spawning_type": "epsilon", "bias_column": 6, "epsilon":0.75, "density": "exitContinuous",
                                      "simulation_type": "pele", "iterations": 1000, "pele_steps": 2, 
                                      "cluster_values": "[1, 2.5]", "cluster_conditions": "[1.1]",
                                      "steric_trials": 500, "overlap_factor": 0.65, "params": pcs.IN_OUT,
                                      "box_radius": 10},
                             "in_out_soft" : {"spawning_type": "independentMetric", "bias_column": 6, "epsilon":0.75, "density": "exitContinuous",
                                      "simulation_type": "pele", "iterations": 1000, "pele_steps": 2, 
                                      "cluster_values": "[1, 2.5]", "cluster_conditions": "[1.1]",
                                      "steric_trials": 500, "overlap_factor": 0.65, "params": pcs.IN_OUT,
                                      "box_radius": 10},
                             "water_exp": {"spawning_type": "independent", "bias_column": 6, "epsilon":0.75, "density": "null",
                                      "simulation_type": "pele", "iterations": 50, "pele_steps": 12, 
                                      "cluster_values": "[1.75, 2.5, 3.5, 5]", "cluster_conditions": "[1.6, 1.2, 1, 0.0]",
                                      "steric_trials": 500, "overlap_factor": 0.65, "params": pcs.WATER_BS,
                                      "box_radius": 10},
                             "water_lig": {"spawning_type": "inverselyProportional", "bias_column": 6,
                                      "epsilon":0.15, "density": "null",
                                      "simulation_type": "pele", "iterations": 50, "pele_steps": 12, 
                                      "cluster_values": "[1.75, 2.5, 3.5, 5]", "cluster_conditions": "[1.6, 1.2, 1, 0.0]",
                                      "steric_trials": 100, "overlap_factor": 0.5, "params": pcs.WATER_LIG,
                                      "box_radius": 10},
                             "bias": {"spawning_type": "epsilon", "bias_column": 5, "epsilon":0.15, "density": "null",
                                      "simulation_type": "pele", "iterations": 50, "pele_steps": 8, 
                                      "cluster_values": "[1.5, 2, 5]", "cluster_conditions": "[0.6, 0.4, 0.0]",
                                      "steric_trials": 250, "overlap_factor": 0.65, "params": pcs.BIAS,
                                      "box_radius": 30},
                             "anm": {"spawning_type": "independent", "bias_column": 5, "epsilon":0.15, "density": "null",
                                      "simulation_type": "pele", "iterations": 50, "pele_steps": 8, 
                                      "cluster_values": "[1.5, 2, 5]", "cluster_conditions": "[0.6, 0.4, 0.0]",
                                      "steric_trials": 250, "overlap_factor": 0.65, "params": pcs.BIAS,
                                      "box_radius": 30}
                    }

             }
        }
        
        software_setings = SOFTWARE_CONSTANTS[args.software]
        if args.software == "msm":
            if args.precision:
                software_setings["files"] = software_setings["files"].get("XP")
            elif args.test:
                software_setings["files"] = software_setings["files"].get("test")
            else:
                software_setings["files"] = software_setings["files"].get("SP")
        if args.software == "glide":
                software_setings["files"] = software_setings["files"].get("glide")
                software_setings["simulation_params"] = software_setings["simulation_params"].get("glide")
        if args.software == "adaptive" :
            if args.full:
                software_setings["files"] = software_setings["files"].get("full")
                software_setings["simulation_params"] = software_setings["simulation_params"].get("full")
            elif args.in_out:
                software_setings["files"] = software_setings["files"].get("in_out")
                software_setings["simulation_params"] = software_setings["simulation_params"].get("in_out")
            elif args.in_out_soft:
                software_setings["files"] = software_setings["files"].get("in_out_soft")
                software_setings["simulation_params"] = software_setings["simulation_params"].get("in_out_soft")
            elif args.induce_fit:
                software_setings["files"] = software_setings["files"].get("induce_fit")
                software_setings["simulation_params"] = software_setings["simulation_params"].get("induce_fit")
            elif args.water_exp:
                software_setings["files"] = software_setings["files"].get("water_exp")
                software_setings["simulation_params"] = software_setings["simulation_params"].get("water_exp")
            elif args.water_lig:
                software_setings["files"] = software_setings["files"].get("water_lig")
                software_setings["simulation_params"] = software_setings["simulation_params"].get("water_lig")
            elif args.bias:
                software_setings["files"] = software_setings["files"].get("bias")
                software_setings["simulation_params"] = software_setings["simulation_params"].get("bias")
            elif args.adaptive and args.pele:
                software_setings["files"] = software_setings["files"].get("adaptive")
                software_setings["simulation_params"] = software_setings["simulation_params"].get("adaptive")
            else:
                #Standard file (user will change the parameters)
                software_setings["files"] = software_setings["files"].get("bias")
                software_setings["simulation_params"] = software_setings["simulation_params"].get("bias")
        
        return software_setings
