"""
This module handles the assignment of optional parameters for FragPELE.
"""


class FragOptionalParameters(object):
    """
    Class to assign the optional parameters for FragPELE.
    """
    def __init__(self, parameters, args):
        """
        Given a Parameters object, it initializes the optional parameters for
        FragPELE.

        .. todo ::
           * We need to unify all classes that prepare the parameters
             for Frag. We need an abstract class to ensure that all of them
             modify correctly the Parameters object

        Parameters
        ----------
        parameters : a Parameters object
            The Parameters object containing the parameters for PELE
        args : a YamlParser object
            The YamlParser object containing the input parameters chosen
            by the user
        """
        from pele_platform.Utilities.Helpers.constraints import \
            alpha_constraints
        from pele_platform.constants import constants

        # Set simulation control file
        parameters.frag_run = args.frag_run

        # Set constraints
        parameters.constraints = alpha_constraints.retrieve_constraints(
            parameters.core, interval=parameters.ca_interval,
            back_constr=parameters.ca_constr,
            ter_constr=parameters.terminal_constr)

        # Set chain
        if args.chain_core:
            parameters.chain_core = args.chain_core
        else:
            parameters.chain_core = parameters.simulation_params.get("chain_core", "L")

        # Set box
        if parameters.box_radius:
            parameters.box = constants.BOX.format(parameters.box_radius,
                                                  parameters.box_center)
        else:
            parameters.box = ""

        # Set ligand
        parameters.gridres = args.gridres
        parameters.plop_path = "PlopRotTemp_S_2017/ligand_prep.py"

        # Set output
        if args.frag_criteria:
            parameters.criteria = args.frag_criteria
        else:
            parameters.criteria = \
                parameters.simulation_params.get("frag_criteria",
                                                 "Binding Energy")

        if args.frag_output_folder:
            parameters.output_folder = args.frag_output_folder
        else:
            parameters.output_folder = \
                parameters.simulation_params.get("frag_output_folder",
                                           "growing_steps")

        if args.frag_cluster_folder:
            parameters.cluster_folder = args.frag_cluster_folder
        else:
            parameters.cluster_folder = \
                parameters.simulation_params.get("frag_cluster_folder",
                                           "clustering_PDBs")

        # Set other parameters
        parameters.distcont = 4
        parameters.threshold = 0.3
        parameters.condition = "min"
        parameters.metricweights = "linear"
        parameters.nclusters = 5
        parameters.min_overlap = 0.5
        parameters.max_overlap = 0.7
        parameters.frag_chain = "L"
        parameters.banned = None
        parameters.limit = None
        parameters.rename = False
        parameters.threshold_clash = 1.7
        parameters.translation_high = 0.05
        parameters.translation_low = 0.02
        parameters.rotation_high = 0.1
        parameters.rotation_low = 0.05
        parameters.explorative = False
        parameters.frag_radius = 10
        parameters.sampling_control = None
        parameters.only_prepare = False
        parameters.only_grow = False
        parameters.no_check = True

        # Interaction restrictions Parameters (only used by adaptive)
        parameters.interaction_restrictions = ""
        parameters.met_interaction_restrictions = ""

        if args.cleanup:
            parameters.cleanup = args.cleanup
        else:
            parameters.cleanup = False
