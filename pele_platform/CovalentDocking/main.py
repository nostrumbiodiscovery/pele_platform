from dataclasses import dataclass
import glob
import os

import pele_platform.Adaptive.simulation as si
import pele_platform.Frag.simulation as fr
import pele_platform.Utilities.Helpers.helpers as hp
import pele_platform.Utilities.Parameters.pele_env as pv
from pele_platform.Utilities.Helpers import bestStructs as bs


@dataclass
class CovalentDocking:
    args: pv.EnviroBuilder

    def run_covalent_docking(self) -> pv.EnviroBuilder:
        self.current_dir = os.getcwd()
        self._set_adaptive_parameters()
        adaptive_simulation = si.run_adaptive(self.args)
        frag_core = self._extract_structures(adaptive_simulation)
        self._set_frag_parameters(adaptive_simulation, frag_core)
        self._write_config_file(adaptive_simulation)
        with hp.cd(adaptive_simulation.pele_dir):
            fragment_growing = fr.FragRunner(self.args).run_simulation()
        
        return adaptive_simulation, fragment_growing

    def _set_adaptive_parameters(self) -> None:

        # Reactive atoms
        self.atom_ligand = self.args.atom_ligand
        self.atom_sidechain = self.args.atom_sidechain  # for example: "A:145"

        # Simulation parameters
        self.args.epsilon = 0.25
        self.args.out_in = True
        self.args.bias_column = '7'  # atom_dist
        self.args.be_column = '7'
        self.args.randomize = True
        self.args.atom_dist = [self.atom_ligand, self.atom_sidechain]

        # Box
        box_center, box_radius = hp.retrieve_box(self.args.system, self.atom_ligand, self.atom_sidechain)
        self.args.box_center = self.args.box_center if self.args.box_center else box_center
        self.args.box_radius = self.args.box_radius if self.args.box_radius else box_radius

    def _extract_structures(self, adaptive_simulation) -> list:

        output_dir = os.path.join(adaptive_simulation.pele_dir, adaptive_simulation.output)
                
        files_out, _, _, _, dist = bs.main(criteria=self.args.be_column, n_structs=1, path=output_dir, topology=None)
        bs_files = [os.path.join(adaptive_simulation.pele_dir, "results", f) for f in files_out]
        
        frag_core = os.path.join(adaptive_simulation.pele_dir, "covalent_docking_input")
        os.system("mkdir {}".format(frag_core))
        
        for file in bs_files:
            os.system("cp {} {}/.".format(file, frag_core))
            hp.remove_residue(os.path.join(frag_core, os.path.basename(file)), adaptive_simulation.residue)

        return frag_core
    
    def _write_config_file(self, adaptive_simulation):
        input_conf = self.args.frag_input
        residue_atom = self.atom_sidechain.split(":")[2]
        ligand_atom = self.atom_ligand.split(":")[2]

        with open(input_conf, "w+") as f:
            f.write("{} {} {}\n".format(adaptive_simulation.ligand_ref, residue_atom, ligand_atom))

        return input_conf

    def _set_frag_parameters(self, adaptive_simulation, frag_core):
        
        # Get rid of adaptive parameters
        self.args.system = None
        self.args.frag_core = glob.glob(os.path.join(frag_core, "*.pdb"))[0]
        self.args.randomize = False
        self.args.out_in = False
        self.args.atom_dist = None

        # Drop atom name from atom string
        self.args.cov_res = self.atom_sidechain.split(":")
        self.args.cov_res = self.args.cov_res[0:2]
        self.args.cov_res = ":".join(self.args.cov_res)
        
        # FragCovalent parameters
        self.args.frag_input = os.path.join(self.current_dir, "input.conf")
        self.be_column = "LocalNonBindingEnergy"

