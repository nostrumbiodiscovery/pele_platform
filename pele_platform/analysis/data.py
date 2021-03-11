"""
This module contains classes and methods to handle data coming from PELE
trajectories.
"""


class DataHandler(object):
    """
    Main class to handle data coming from PELE trajectories.
    """

    def __init__(self, parameters):
        """
        It initializes a DataHandler object.
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

    def get_reports_dataframe(self):
        """
        It returns the data stored in PELE reports as a pandas dataframe.

        Returns
        -------
        dataframe : a pandas.DataFrame object
            The dataframe containing the information from PELE reports
        """
        import os
        import glob
        from AdaptivePELE.utilities import utilities
        import pandas as pd
        from pele_platform.Utilities.Helpers import get_suffix

        # Initialize primary variables
        sim_path = os.path.join(self.parameters.pele_dir,
                                self.parameters.output)
        epoch_dirs = glob.glob(os.path.join(sim_path, '[0-9]*'))
        report_prefix = self.parameters.report_name
        trajectory_prefix = \
            str(os.path.splitext(self.parameters.traj_name)) + '_'
        trajectory_format = \
            str(os.path.splitext(self.parameters.traj_name)[-1])

        # Filter out non digit folders
        epochs = [os.path.basename(path) for path in epoch_dirs
                  if os.path.basename(path).isdigit()]

        dataframe_lists = []
        for adaptive_epoch in sorted(epochs, key=int):
            folder = os.path.join(sim_path, str(adaptive_epoch))
            report_dirs = glob.glob(os.path.join(folder,
                                                 report_prefix + '_[0-9]*'))

            report_ids = [get_suffix(path) for path in report_dirs
                          if get_suffix(path).isdigit()]
            report_list = [os.path.join(folder, report_prefix + '_' + i)
                           for i in sorted(report_ids, key=int)]

            for i, report in enumerate(report_list, start=1):
                pandas_df = pd.read_csv(report, sep="    ", engine="python",
                                        index_col=False, header=0)
                pandas_df["epoch"] = adaptive_epoch
                pandas_df["trajectory"] = \
                    os.path.join(sim_path, adaptive_epoch,
                                 trajectory_prefix + '_' + str(i) +
                                 trajectory_format)
                dataframe_lists.append(pandas_df)

        dataframe = pd.concat(dataframe_lists, ignore_index=True)

        return dataframe
