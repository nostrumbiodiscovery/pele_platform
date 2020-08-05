import pele_platform.Adaptive.simulation as si
import pele_platform.Utilities.Helpers.helpers as hp
import pele_platform.Utilities.Parameters.pele_env as pv


class CovalentDocking:
    args: pv.EnviroBuilder

    def run_covalent_docking(self) -> pv.EnviroBuilder:
        self._set_parameters()
        simulation_parameters = si.run_adaptive(self.args)

        return simulation_parameters

    def _set_parameters(self) -> None:

        # Simulation parameters
        self.epsilon = 0.25
        self.out_in = True
        self.bias_column = 7  # atom_dist
        self.args.randomize = True

        # Reactive atoms
        self.atom_ligand = self.args.atom_ligand
        self.atom_sidechain = self.args.atom_sidechain

        # Box
        box_center, box_radius = hp.retrieve_box(self.args.system, self.atom_ligand, self.atom_sidechain)
        self.args.box_center = self.args.box_center if self.args.box_center else box_center
        self.args.box_radius = self.args.box_radius if self.args.box_radius else box_radius
