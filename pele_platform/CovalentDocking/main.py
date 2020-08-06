from dataclasses import dataclass
import pele_platform.Adaptive.simulation as si
import pele_platform.Utilities.Helpers.helpers as hp
import pele_platform.Utilities.Parameters.pele_env as pv
from pele_platform.Utilities.Helpers import bestStructs as bs
import os

@dataclass
class CovalentDocking:
    args: pv.EnviroBuilder

    def run_covalent_docking(self) -> pv.EnviroBuilder:
        self._set_adaptive_parameters()
        adaptive_simulation = si.run_adaptive(self.args)
        frag_input = self._extract_structures(adaptive_simulation)
        self._set_frag_parameters(adaptive_simulation, frag_input)
        print("adaptive pele dir", adaptive_simulation.pele_dir) 
        with hp.cd(adaptive_simulation.pele_dir):
            fragment_growing = si.run_adaptive(self.args)
        
        return adaptive_simulation, fragment_growing

    def _set_adaptive_parameters(self) -> None:

        # Reactive atoms
        self.atom_ligand = self.args.atom_ligand
        self.atom_sidechain = self.args.atom_sidechain

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
        
        frag_input = os.path.join(adaptive_simulation.pele_dir, "covalent_docking_input")
        os.system("mkdir {}".format(frag_input))
        
        for file in bs_files:
            os.system("cp {} {}/.".format(file, frag_input))
            hp.remove_residue(os.path.join(frag_input, os.path.basename(file)), adaptive_simulation.residue)

        return frag_input
    
    def _set_frag_parameters(self, adaptive_simulation, frag_input):
        
        # Get rid of adaptive parameters
        self.args.system = os.path.join(frag_input, "*.pdb")
        self.args.resname = None
        self.args.randomize = False
        self.args.out_in = False
        self.args.atom_dist = None

        # Set up frag
        
