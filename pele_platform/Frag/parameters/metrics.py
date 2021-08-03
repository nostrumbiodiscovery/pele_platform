"""
This module handles the assignment of metrics parameters for FragPELE.
"""


class FragMetrics(object):
    """
    Base class to assign the metrics parameters for FragPELE.
    """

    def __init__(self, parameters, args):
        """
        Given a Parameters object, it initializes the metrics parameters for
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
        from pele_platform.Adaptive.metrics import MetricBuilder

        metric_builder = MetricBuilder()

        if args.atom_dist:
            parameters.metrics = \
                metric_builder.distance_to_atom_json(parameters.system,
                                                     args.atom_dist)
        else:
            parameters.metrics = ""

        if args.native:
            parameters.native = metric_builder.rmsd_to_json(args.native,
                                                            parameters.chain)
        else:
            parameters.native = ""

        parameters.local_nonbonding_energy = ""
