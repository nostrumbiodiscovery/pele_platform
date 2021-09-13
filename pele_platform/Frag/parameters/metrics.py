"""
This module handles the assignment of metrics parameters for FragPELE.
"""
from pele_platform.context import context


class FragMetrics(object):
    """
    Base class to assign the metrics parameters for FragPELE.
    """

    def __init__(self):
        """
        Given a Parameters object, it initializes the metrics parameters for
        FragPELE.

        .. todo ::
           * We need to unify all classes that prepare the parameters
             for Frag. We need an abstract class to ensure that all of them
             modify correctly the Parameters object
        """
        from pele_platform.Adaptive.metrics import MetricBuilder

        metric_builder = MetricBuilder()

        if context.yaml_parser.atom_dist:
            context.parameters.metrics = \
                metric_builder.distance_to_atom_json(context.parameters.system,
                                                     context.yaml_parser.atom_dist)
        else:
            context.parameters.metrics = ""

        if context.yaml_parser.native:
            context.parameters.native = metric_builder.rmsd_to_json(context.yaml_parser.native,
                                                                    context.parameters.chain)
        else:
            context.parameters.native = ""

        context.parameters.local_nonbonding_energy = ""
