import os
import pele_platform.constants.constants as cs
import pele_platform.constants.pele_params as pcs

"""
Description of the file: This file specifies what files to use
for each job.

"""


def retrieve_software_settings(args, pele_dir):

    SOFTWARE_CONSTANTS = {
        "adap_ex_output": None,
        "ad_ex_temp": os.path.join(pele_dir, "adaptive.conf"),
        "pele_exit_temp": os.path.join(pele_dir, "pele.conf"),
        "folders": [
            "",
            "DataLocal/Templates/OPLS2005/HeteroAtoms/",
            "DataLocal/Templates/AMBER99sb/HeteroAtoms/",
            "DataLocal/Templates/AMBER99sbBSC0/HeteroAtoms/",
            "DataLocal/LigandRotamerLibs",
            "DataLocal/OBC",
        ],
        "file_names": ["adaptive.conf", "pele.conf"],
        "files": [
            os.path.join(cs.DIR, "Templates/template_adaptive.conf"),
            os.path.join(cs.DIR, "Templates/pele_template.conf"),
        ],
        "simulation_params": {
            "out_in": {
                "spawning_type": "inverselyProportional",
                "density": "continuous",
                "simulation_type": "pele",
                "iterations": 100,
                "pele_steps": 8,
                "cluster_values": "[2, 5, 7]",
                "cluster_conditions": "[1, 0.6, 0.0]",
                "steric_trials": 250,
                "overlap_factor": 0.65,
                "params": pcs.OUT_IN,
                "box_radius": 30,
            },
            "gpcr_orth": {
                "spawning_type": "epsilon",
                "density": "null",
                "simulation_type": "pele",
                "iterations": 50,
                "pele_steps": 8,
                "cluster_values": "[1.75, 2.5, 4]",
                "cluster_conditions": "[0.7, 0.4, 0.0]",
                "steric_trials": 100,
                "overlap_factor": 0.65,
                "epsilon": 0.25,
                "bias_column": 6,
                "params": pcs.GPCR_ORTH,
                "box_radius": 15,
                "ca_constr": 5.0,
                "terminal_constr": 5.0,
                "ca_interval": 5,
            },
            "global": {
                "spawning_type": "inverselyProportional",
                "bias_column": 5,
                "epsilon": 0.25,
                "density": "continuous",
                "simulation_type": "pele",
                "iterations": 100,
                "pele_steps": 8,
                "cluster_values": "[2.0, 5, 7]",
                "cluster_conditions": "[1, 0.6, 0.0]",
                "steric_trials": 200,
                "overlap_factor": 0.65,
                "params": pcs.GLOBAL,
                "box_radius": None,
            },
            "induced_fit_exhaustive": {
                "spawning_type": "independent",
                "bias_column": 5,
                "epsilon": 0.25,
                "density": "null",
                "simulation_type": "pele",
                "iterations": 1,
                "pele_steps": 1000,
                "cluster_values": "[2.0, 5, 7]",
                "cluster_conditions": "[2.0, 5, 7]",
                "steric_trials": 500,
                "overlap_factor": 0.65,
                "params": pcs.INDUCED_FIT,
                "box_radius": 6,
            },
            "induced_fit_fast": {
                "spawning_type": "inverselyProportional",
                "bias_column": 5,
                "epsilon": 0.25,
                "density": "null",
                "simulation_type": "pele",
                "iterations": 30,
                "pele_steps": 12,
                "cluster_values": "[2.0, 5, 7]",
                "cluster_conditions": "[2.0, 5, 7]",
                "steric_trials": 500,
                "overlap_factor": 0.65,
                "params": pcs.INDUCED_FIT,
                "box_radius": 6,
            },
            "in_out": {
                "spawning_type": "epsilon",
                "bias_column": 6,
                "epsilon": 0.75,
                "density": "exitContinuous",
                "simulation_type": "pele",
                "iterations": 1000,
                "pele_steps": 2,
                "cluster_values": "[1, 2.5]",
                "cluster_conditions": "[1.1]",
                "steric_trials": 500,
                "overlap_factor": 0.65,
                "params": pcs.IN_OUT,
                "box_radius": 10,
            },
            "in_out_soft": {
                "spawning_type": "independentMetric",
                "bias_column": 6,
                "epsilon": 0.75,
                "density": "exitContinuous",
                "simulation_type": "pele",
                "iterations": 1000,
                "pele_steps": 2,
                "cluster_values": "[1, 2.5]",
                "cluster_conditions": "[1.1]",
                "steric_trials": 500,
                "overlap_factor": 0.65,
                "params": pcs.IN_OUT,
                "box_radius": 10,
                "spawning_condition": "max",
            },
            "rescoring": {
                "spawning_type": "independent",
                "bias_column": 5,
                "epsilon": 0.25,
                "density": "null",
                "simulation_type": "pele",
                "iterations": 20,
                "pele_steps": 12,
                "cluster_values": "[1.75, 2.5, 4, 6]",
                "cluster_conditions": "[1, 0.6, 0.4, 0.0]",
                "steric_trials": 500,
                "overlap_factor": 0.65,
                "params": pcs.RESCORING,
                "box_radius": 6,
                "anm_freq": 6,
                "sidechain_freq": 3,
                "min_freq": 1,
                "temperature": 1000,
                "anm_displacement": 0.5,
                "anm_modes_change": 3,
                "ca_constr": 2.5,
                "terminal_constr": 5.0,
                "ca_interval": 8,
            },
            "anm": {
                "spawning_type": "independent",
                "bias_column": 5,
                "epsilon": 0.15,
                "density": "null",
                "simulation_type": "pele",
                "iterations": 50,
                "pele_steps": 8,
                "cluster_values": "[1.5, 2, 5]",
                "cluster_conditions": "[0.6, 0.4, 0.0]",
                "steric_trials": 250,
                "overlap_factor": 0.65,
                "params": pcs.OUT_IN,
                "box_radius": 30,
            },
        },
    }

    software_setings = SOFTWARE_CONSTANTS
    if args.full:
        type_simulation = "global"
    elif args.in_out:
        type_simulation = "in_out"
    elif args.in_out_soft:
        type_simulation = "in_out_soft"
    elif args.induced_fit_exhaustive:
        type_simulation = "induced_fit_exhaustive"
    elif args.induced_fit_fast:
        type_simulation = "induced_fit_fast"
    elif args.rescoring:
        type_simulation = "rescoring"
    elif args.out_in:
        type_simulation = "out_in"
    elif args.gpcr_orth:
        type_simulation = "gpcr_orth"
    else:
        # Standard file (user will change the parameters)
        type_simulation = "induced_fit_fast"

    software_setings["files"] = software_setings["files"]
    software_setings["simulation_params"] = software_setings["simulation_params"].get(
        type_simulation
    )
    return software_setings
