import os
import pele_platform.constants as cs

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
                                 }
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
                                         os.path.join(cs.DIR, "Templates/adaptive_glide.conf"),
                                         os.path.join(cs.DIR, "Templates/pele_glide.conf")],
                              "test" : [os.path.join(cs.DIR, "Templates/glide.in"),
                                        os.path.join(cs.DIR, "Templates/adaptive_glide_test.conf"),
                                        os.path.join(cs.DIR, "Templates/pele_glide.conf")],
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
                              "out_in" : [ os.path.join(cs.DIR, "Templates/out_in_adaptive.conf"), 
                                           os.path.join(cs.DIR, "Templates/out_in_pele.conf") ],
                              "full" : [ os.path.join(cs.DIR, "Templates/global_adaptive.conf"), 
                                           os.path.join(cs.DIR, "Templates/global_pele.conf") ],
                              "induce_fit" : [ os.path.join(cs.DIR, "Templates/induce_fit_adaptive.conf"),
                                               os.path.join(cs.DIR, "Templates/induce_fit_pele.conf") ],
                              "in_out" : [ os.path.join(cs.DIR, "Templates/in_out_adaptive.conf"),
                                           os.path.join(cs.DIR, "Templates/in_out_pele.conf") ],
                              "in_out_soft" : [ os.path.join(cs.DIR, "Templates/in_out_soft_adaptive.conf"),
                                           os.path.join(cs.DIR, "Templates/in_out_pele.conf") ],
                              "water_exp":  [ os.path.join(cs.DIR, "Templates/adaptive_water.conf"),
                                           os.path.join(cs.DIR, "Templates/pele_water.conf") ],
                              "water_lig":  [ os.path.join(cs.DIR, "Templates/adaptive_water_ligand.conf"),
                                           os.path.join(cs.DIR, "Templates/pele_water_ligand.conf") ],
                              "bias":  [ os.path.join(cs.DIR, "Templates/biased_adaptive.conf"),
                                           os.path.join(cs.DIR, "Templates/biased_pele.conf") ]
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
            if args.test:
                software_setings["files"] = software_setings["files"].get("test")
            else:
                software_setings["files"] = software_setings["files"].get("glide")
        if args.software == "adaptive" :
            if args.out_in:
                software_setings["files"] = software_setings["files"].get("out_in")
            elif args.full:
                software_setings["files"] = software_setings["files"].get("full")
            elif args.in_out:
                software_setings["files"] = software_setings["files"].get("in_out")
            elif args.in_out_soft:
                software_setings["files"] = software_setings["files"].get("in_out_soft")
            elif args.induce_fit:
                software_setings["files"] = software_setings["files"].get("induce_fit")
            elif args.water_exp:
                software_setings["files"] = software_setings["files"].get("water_exp")
            elif args.water_lig:
                software_setings["files"] = software_setings["files"].get("water_lig")
            elif args.bias:
                software_setings["files"] = software_setings["files"].get("bias")
            else:
                software_setings["files"] = software_setings["files"].get("adaptive")
        
        return software_setings
