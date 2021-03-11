"""
This module manages the analysis toolkit of the platform.
"""


__all__ = ["Analysis"]


class Analysis(object):
    """
    General class to manage all analysis operations.
    """

    def __init__(self, parameters):
        """
        It initializes an Analysis instance which it depends on
        the general Parameters class of the PELE Platform.

        Parameters
        ----------
        parameters : a Parameters object
            The Parameters object containing the parameters that belong
            to the simulation
        """
        self._parameters = parameters

    @property
    def parameters(self):
        """
        It returns the Parameters object to analyze.

        Returns
        -------
        parameters : a Parameters object
            The Parameters object containing the parameters that belong
            to the simulation
        """
        return self._parameters
