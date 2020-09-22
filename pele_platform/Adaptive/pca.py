from dataclasses import dataclass
import glob
import os
import pele_platform.constants.constants as cs
import pele_platform.Utilities.Helpers.helpers as helpers
import pele_platform.Utilities.Helpers.calculatePCA4PELE as pc


@dataclass
class PCA:
    """
    Base class to perform pca analysis
    on user's receptor traj and generate
    the json input's code. Two input options:

    1) pca_traj = trajectory*.pdb
    2) pca_traj = [trajectory1.pdb, trajectory2.pdb, ...]
    """

    pca_traj: list
    pele_dir: str

    def generate(self, env) -> str:
        # Calculate pca and retrieve json format
        pdbs = self._retrieve_trajectories()
        pca_json = self._pca_to_json(pdbs, logger=env)
        return pca_json

    def _retrieve_trajectories(self) -> list:
        # Get pdb inputs to build pca
        if isinstance(self.pca_traj, str):
            pdbs = glob.glob(self.pca_traj)
        elif isinstance(self.pca_traj, list):
            pdbs = self.pca_traj
        return pdbs

    def _pca_to_json(self, pdbs, logger=None) -> str:
        # calculate pca over pdb and get json
        pdbs_full_path = [os.path.abspath(pdb) for pdb in pdbs]
        with helpers.cd(self.pele_dir):
            pca = pc.main(pdbs_full_path, logger=logger)
        self.pca = cs.PCA.format(pca)
        return self.pca
