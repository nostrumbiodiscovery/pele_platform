import os
import pele_platform.constants.constants as cs


class SimulationPaths:
    def __init__(self):
        self.pele_dir = self.pele_dir
        self.inputs_dir = os.path.join(self.pele_dir, "input")
        self.output = os.path.join(self.pele_dir, "output")

        if self.no_ppp:
            self.system_fix = os.path.join(
                self.inputs_dir, os.path.basename(self.system)
            )
        else:
            self.system_fix = os.path.join(
                self.inputs_dir,
                "{}_processed.pdb".format(
                    os.path.splitext(os.path.basename(self.system))[0]
                ),
            )

        self.adap_ex_input = os.path.basename(self.system_fix)
        self.pele_temp = os.path.join(self.pele_dir, "pele.conf")
        self.adap_l_input = "{}/initial_*"
        self.adap_l_output = os.path.join(self.pele_dir, "output_pele")
        self.ad_l_temp = os.path.join(self.pele_dir, "adaptive_long.conf")
        self.ligand_ref = os.path.join(self.inputs_dir, "ligand.pdb")
        self.rotamers_folder = os.path.join(
            self.pele_dir, "DataLocal/LigandRotamerLibs/"
        )
        self.template_folder = os.path.join(
            self.pele_dir, f"DataLocal/Templates/{self.forcefield}/HeteroAtoms/"
        )
        self.receptor = os.path.join(self.inputs_dir, "receptor.pdb")
        self.topology = (
            None
            if self.pdb
            else os.path.join(self.pele_dir, self.output, "topologies/topology_0.pdb")
        )
        self.obc_tmp = os.path.join(cs.DIR, "Templates/solventParamsHCTOBC.txt")
        self.obc_file = os.path.join(
            self.pele_dir, "DataLocal/OBC/solventParamsHCTOBC.txt"
        )
        self.box_temp = os.path.join(self.pele_dir, "box.pdb")

        # shit from features/adaptive SOFTWARE CONSTANTS
        self.adap_ex_output = None
        self.ad_ex_temp = os.path.join(self.pele_dir, "adaptive.conf")
        self.pele_exit_temp = os.path.join(self.pele_dir, "pele.conf")
        self.folders = [
            "",
            "DataLocal/Templates/OPLS2005/HeteroAtoms/",
            "DataLocal/Templates/AMBER99sb/HeteroAtoms/",
            "DataLocal/Templates/AMBER99sbBSC0/HeteroAtoms/",
            "DataLocal/LigandRotamerLibs",
            "DataLocal/OBC",
        ]
        self.file_names = [
            "adaptive.conf",
            "pele.conf",
        ]
        self.files = [
                         os.path.join(cs.DIR, "Templates/template_adaptive.conf"),
                         os.path.join(cs.DIR, "Templates/pele_template.conf"),
                     ]
