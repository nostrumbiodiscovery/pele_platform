"""
This module handles the assignment of simulation parameters for FragPELE.
"""
from pele_platform.constants import constants as cs
from pele_platform.context import context


class FragSimulationParameters(object):
    """
    Base class to assign the simulation parameters for FragPELE.
    """

    def __init__(self):
        """
        Given a Parameters object, it initializes the simulation parameters for
        FragPELE.

        .. todo ::
           * We need to unify all classes that prepare the parameters
             for Frag. We need an abstract class to ensure that all of them
             modify correctly the Parameters object
        """
        import os

        from pele_platform.constants import constants

        context.parameters.system = context.yaml_parser.frag_core

        # Set number of growing steps
        if context.yaml_parser.growing_steps:
            context.parameters.gr_steps = context.yaml_parser.growing_steps
        else:
            context.parameters.gr_steps = \
                context.parameters.simulation_params.get("growing_steps", 6)

        # Set number of FragPELE steps
        if context.yaml_parser.frag_steps:
            context.parameters.frag_steps = context.yaml_parser.frag_steps
        else:
            context.parameters.frag_steps = \
                context.parameters.simulation_params.get("steps_in_gs", 3)

        # Set number of equilibration steps
        if context.yaml_parser.frag_eq_steps:
            context.parameters.frag_eq_steps = context.yaml_parser.frag_eq_steps
        else:
            context.parameters.frag_eq_steps = \
                context.parameters.simulation_params.get("sampling_steps", 20)

        # Set path to control file
        context.parameters.control_file = os.path.join(constants.DIR,
                                                       "Templates/pele_template.conf")

        # Set FragPELE protocol
        if context.yaml_parser.protocol:
            context.parameters.protocol = context.yaml_parser.protocol
        else:
            context.parameters.simulation_params.get("protocol", "")

        # Set path to topology file
        if context.parameters.pdb:
            context.parameters.topology = None
        else:
            context.parameters.topology = os.path.join("output_pele", "topology.pdb")

        # Set FragPELE restart option
        if context.yaml_parser.frag_restart:
            context.parameters.frag_restart = "-rst"
        else:
            context.parameters.frag_restart = ""

        # Set trajectory name
        context.parameters.frag_traj_name = "trajectory"

        # Set use srun parameter for FragPELE
        if context.parameters.usesrun == "true":
            context.parameters.usesrun = True
        else:
            context.parameters.usesrun = False

        context.parameters.spython = cs.SCHRODINGER
