"""
This module handles the assignment of optional parameters for FragPELE.
"""
from pele_platform.context import context


class FragOptionalParameters(object):
    """
    Class to assign the optional parameters for FragPELE.
    """

    def __init__(self):
        """
        Given a Parameters object, it initializes the optional parameters for
        FragPELE.

        .. todo ::
           * We need to unify all classes that prepare the parameters
             for Frag. We need an abstract class to ensure that all of them
             modify correctly the Parameters object
        """
        from pele_platform.Utilities.Helpers.constraints import \
            alpha_constraints
        from pele_platform.constants import constants

        # Set simulation control file
        context.parameters.frag_run = context.yaml_parser.frag_run

        # Set constraints
        context.parameters.constraints = alpha_constraints.retrieve_constraints(
            context.parameters.core, interval=context.parameters.ca_interval,
            back_constr=context.parameters.ca_constr,
            ter_constr=context.parameters.terminal_constr)

        # Set chain
        if context.yaml_parser.chain_core:
            context.parameters.chain_core = context.yaml_parser.chain_core
        else:
            context.parameters.chain_core = context.parameters.simulation_params.get("chain_core", "L")

        # Set box
        if context.parameters.box_radius:
            context.parameters.box = constants.BOX.format(context.parameters.box_radius,
                                                          context.parameters.box_center)
        else:
            context.parameters.box = ""

        # Set ligand
        context.parameters.gridres = context.yaml_parser.gridres
        context.parameters.plop_path = "PlopRotTemp_S_2017/ligand_prep.py"

        # Set output
        if context.yaml_parser.frag_criteria:
            context.parameters.criteria = context.yaml_parser.frag_criteria
        else:
            context.parameters.criteria = \
                context.parameters.simulation_params.get("frag_criteria",
                                                         "Binding Energy")

        if context.yaml_parser.frag_output_folder:
            context.parameters.output_folder = context.yaml_parser.frag_output_folder
        else:
            context.parameters.output_folder = \
                context.parameters.simulation_params.get("frag_output_folder",
                                                         "growing_steps")

        if context.yaml_parser.frag_cluster_folder:
            context.parameters.cluster_folder = context.yaml_parser.frag_cluster_folder
        else:
            context.parameters.cluster_folder = \
                context.parameters.simulation_params.get("frag_cluster_folder",
                                                         "clustering_PDBs")

        # Set other parameters
        context.parameters.distcont = 4
        context.parameters.threshold = 0.3
        context.parameters.condition = "min"
        context.parameters.metricweights = "linear"
        context.parameters.nclusters = 5
        context.parameters.min_overlap = 0.5
        context.parameters.max_overlap = 0.7
        context.parameters.frag_chain = "L"
        context.parameters.banned = None
        context.parameters.limit = None
        context.parameters.rename = False
        context.parameters.threshold_clash = 1.7
        context.parameters.translation_high = 0.05
        context.parameters.translation_low = 0.02
        context.parameters.rotation_high = 0.1
        context.parameters.rotation_low = 0.05
        context.parameters.explorative = False
        context.parameters.frag_radius = 10
        context.parameters.sampling_control = None
        context.parameters.only_prepare = False
        context.parameters.only_grow = False
        context.parameters.no_check = True

        # Interaction restrictions Parameters (only used by adaptive)
        context.parameters.interaction_restrictions = ""
        context.parameters.met_interaction_restrictions = ""

        if context.yaml_parser.cleanup:
            context.parameters.cleanup = context.yaml_parser.cleanup
        else:
            context.parameters.cleanup = False
