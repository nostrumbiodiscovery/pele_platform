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
        self._set_parameters()
        simulation_parameters = si.run_adaptive(self.args)
        frag_input = self._extract_structures(simulation_parameters)
        
        return simulation_parameters

    def _set_parameters(self) -> None:

        # Reactive atoms
        self.atom_ligand = self.args.atom_ligand
        self.atom_sidechain = self.args.atom_sidechain

        # Simulation parameters
        self.epsilon = 0.25
        self.out_in = True
        self.args.bias_column = '7'  # atom_dist
        self.args.be_column = '7'
        self.args.randomize = True
        self.args.atom_dist = [self.atom_ligand, self.atom_sidechain]

        # Box
        box_center, box_radius = hp.retrieve_box(self.args.system, self.atom_ligand, self.atom_sidechain)
        self.args.box_center = self.args.box_center if self.args.box_center else box_center
        self.args.box_radius = self.args.box_radius if self.args.box_radius else box_radius

    def _extract_structures(self, simulation_parameters) -> list:

        max_dist = 2.0 # might need to change that
        output_dir = os.path.join(simulation_parameters.pele_dir, simulation_parameters.output)
        self.n_structs = 5 if simulation_parameters.test else 100
        
        files_out, _, _, _, dist = bs.main(criteria=self.args.be_column, n_structs=self.n_structs, path=output_dir, topology=None)
        frag_input = [f for f, d in zip(files_out, dist) if d <= max_dist]

        return frag_input
