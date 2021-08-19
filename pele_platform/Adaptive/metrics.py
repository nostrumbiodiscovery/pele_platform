import pele_platform.Utilities.Helpers.helpers as hp
import pele_platform.constants.constants as cs
import os


class MetricBuilder:
    def __init__(self):
        self.distance = None
        self.rmsd = None
        self.local_nonbonding_energy = None

    def distance_to_atom_json(self, pdb, atom_dist):
        distances = []

        for i in range(1, len(atom_dist), 2):
            atom1, atom1_tag = hp.parse_atom_dist(atom_dist[i - 1], pdb)
            atom2, atom2_tag = hp.parse_atom_dist(atom_dist[i], pdb)
            distances.append(
                cs.DISTANCE_ATOMS.format(
                    atom1=atom1,
                    atom1_tag=atom1_tag,
                    atom2=atom2,
                    atom2_tag=atom2_tag,
                    metric_name=i / 2,
                )
            )
        self.distance = self._distance_to_json(distances)
        return self.distance

    def _distance_to_json(self, distances):
        return "\n".join(distances)

    def rmsd_to_json(self, pdb_reference, chain_rmsd):
        self.rmsd = cs.NATIVE.format(os.path.abspath(pdb_reference), chain_rmsd)
        return self.rmsd

    def local_nonbonding_energy_json(self, covalent_residue, radius):
        if covalent_residue:
            self.local_nonbonding_energy = cs.LOCAL_NONBONDING_ENERGY.format(
                covalent_residue, radius
            )
        else:
            self.local_nonbonding_energy = ""
        return self.local_nonbonding_energy
