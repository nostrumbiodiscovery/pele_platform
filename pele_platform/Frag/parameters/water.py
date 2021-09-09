"""
This module handles the assignment of water parameters for FragPELE.
"""
from pele_platform.context import context


class FragWaterParams(object):
    """
    Class to assign the water parameters for FragPELE.
    """

    def __init__(self):
        """
        Given a Parameters object, it initializes the water parameters for
        FragPELE.

        .. todo ::
           * We need to unify all classes that prepare the parameters
             for Frag. We need an abstract class to ensure that all of them
             modify correctly the Parameters object
        """
        # We first add a new attribute 'pdb' to the Parameters object
        # pointing to the frag_core path
        context.parameters.pdb = context.yaml_parser.frag_core

        # We then initialize the water parameters
        self._initialize_waters()

    @staticmethod
    def _initialize_waters():
        """
        It initializes the water context.parameters.
        """
        from pele_platform.Utilities.Helpers import water

        context.parameters.water_object = water.WaterIncluder(
            [context.parameters.pdb], context.parameters.n_waters,
            user_waters=context.parameters.waters,
            ligand_perturbation_params=context.parameters.parameters,
            water_center=context.parameters.water_center,
            water_radius=context.parameters.water_radius,
            allow_empty_selectors=context.parameters.allow_empty_selectors,
            water_temp=context.parameters.water_temp,
            water_trials=context.parameters.water_trials,
            water_overlap=context.parameters.water_overlap,
            water_constr=context.parameters.water_constr,
            test=context.parameters.test,
            water_freq=context.parameters.water_freq,
            ligand_residue=context.parameters.residue)

        context.parameters.water_object.run()

        # Save water IDs as main parameters attribute for analysis
        context.parameters.water_ids_to_track = context.parameters.water_object.water_ids_to_track
