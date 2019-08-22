import pele_platform.Utilities.Helpers.helpers as hp
import pele_platform.constants.constants as cs


class Metrics_Builder():


    def __init__(self, pdb):
        self.pdb = pdb
        self.metrics = []

    def distance_to_atom(self, atom_dist):
        for i in range(1, len(atom_dist), 2):
            atom1 = hp.retrieve_atom_info(atom_dist[i-1], self.pdb)
            atom2 = hp.retrieve_atom_info(atom_dist[i], self.pdb)
            self.metrics.append(cs.DISTANCE_ATOMS.format(atom1, atom2, i/2))

    def get_metrics(self):
        if self.metrics:
            return self.metrics
        else:
            return ["",]
