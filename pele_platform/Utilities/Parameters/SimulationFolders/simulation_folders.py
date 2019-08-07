import random
import os
import pele_platform.constants as cs
from pele_platform.Utilities.Parameters.SimulationFolders.MSMFolders import msm_folders
from pele_platform.Utilities.Parameters.SimulationFolders.GlideFolders import glide_folders
import pele_platform.Utilities.Helpers.helpers as hp


class SimulationPaths(msm_folders.MSMPaths, glide_folders.GlidePaths):


    def __init__(self, args):
        self.working_folder_paths(args)
        self.ligand_paths(args)
        self.complex_paths(args)
        self.solvent_paths(args)
        self.box_paths(args)

        if self.software == "msm":
            msm_folders.MSMPaths.__init__(self, args)
        elif self.software == "glide":
            glide_folders.GlidePaths.__init__(self, args)


    def working_folder_paths(self, args):
        pele_dir = os.path.abspath("{}_Pele".format(self.residue))

        if not self.folder:
            self.pele_dir = hp.is_repited(pele_dir) if self.restart in cs.FIRST_RESTART else hp.is_last(pele_dir)
        else:
            self.pele_dir = os.path.abspath(self.folder)

        if self.mae_lig:
            self.system_fix = os.path.join(self.pele_dir, "{}_complex_processed.pdb".format(os.path.splitext(os.path.basename(self.system))[0]))
        else:
            self.system_fix = os.path.join(self.pele_dir, "{}_processed.pdb".format(os.path.splitext(os.path.basename(self.system))[0]))

        self.adap_ex_input = os.path.join(self.pele_dir, os.path.basename(self.system_fix))
        self.pele_temp = os.path.join(self.pele_dir, "pele.conf")
        self.adap_l_input = "{}/initial_*"
        self.adap_l_output = os.path.join(self.pele_dir, "output_pele")
        self.ad_l_temp = os.path.join(self.pele_dir, "adaptive_long.conf")

    def ligand_paths(self, args):
        self.ligand_ref = os.path.join(self.pele_dir, "ligand.pdb")
        self.rotamers_folder = os.path.join(self.pele_dir, "DataLocal/LigandRotamerLibs/")
        self.template_folder = os.path.join(self.pele_dir, "DataLocal/Templates/{}/HeteroAtoms/".format(self.forcefield))

    def complex_paths(self, args):
        self.receptor = os.path.join(self.pele_dir, "receptor.pdb")
        self.topology = None if self.pdb else os.path.join("output_pele", "topology.pdb")

    def solvent_paths(self, args):
        self.obc_tmp = os.path.join(cs.DIR, "Templates/solventParamsHCTOBC.txt")
        self.obc_file = os.path.join(self.pele_dir, "DataLocal/OBC/solventParamsHCTOBC.txt")

    def box_paths(self, args):
        self.box_temp = os.path.join(self.pele_dir, "box.pdb")
