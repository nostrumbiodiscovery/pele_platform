"""
This module handles the assignment of water parameters for FragPELE.
"""


class FragWaterParams(object):
    """
    Class to assign the water parameters for FragPELE.
    """

    def __init__(self, parameters, args):
        """
        Given a Parameters object, it initializes the water parameters for
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
        # We first add a new attribute 'pdb' to the Parameters object
        # pointing to the frag_core path
        parameters.pdb = args.frag_core

        # We then initialize the water parameters
        self._initialize_waters(parameters)

    @staticmethod
    def _initialize_waters(parameters):
        """
        It initializes the water parameters.

        Parameters
        ----------
        parameters : a Parameters object
            The Parameters object containing the parameters for PELE
        """
        from pele_platform.Utilities.Helpers import water

        parameters.water_object = water.WaterIncluder(
            [parameters.pdb], parameters.n_waters,
            user_waters=parameters.waters,
            ligand_perturbation_params=parameters.parameters,
            water_center=parameters.water_center,
            water_radius=parameters.water_radius,
            allow_empty_selectors=parameters.allow_empty_selectors,
            water_temp=parameters.water_temp,
            water_trials=parameters.water_trials,
            water_overlap=parameters.water_overlap,
            water_constr=parameters.water_constr,
            test=parameters.test,
            water_freq=parameters.water_freq,
            ligand_residue=parameters.residue)

        parameters.water_object.run()

        # Save water IDs as main parameters attribute for analysis
        parameters.water_ids_to_track = parameters.water_object.water_ids_to_track
