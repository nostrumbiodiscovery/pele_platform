import pele_platform.Utilities.Helpers.helpers as hp
import pele_platform.constants.constants as cs
import os


class StructuralConstraintsBuilder():

    def __init__(self):
        self.metrics = []
        self.conditions = []

    def parse_structural_metrics (self, pdb, constraints_conf):
        for i in range (0, len(constraints_conf)):
            actual = constraints_conf[i]
            name = ""
            if 'name' in actual:
                name = actual['name']
            if 'condition' in actual:
                self.conditions.append(actual['condition'])
            if 'type' in actual:
                self._add_metric (pdb, actual['type'], actual['values'], name)

    def conditions_to_json (self):
        return cs.STRUCTURAL_CONSTRAINTS.format('",\n\t"'.join(self.conditions))

    def metrics_to_json (self):
        return "\n".join(self.metrics)

    def _add_metric (self, pdb, type, values, name):
        if (type == "atom_dist"):
                atom1 = hp.retrieve_atom_info(values[0], pdb)
                atom2 = hp.retrieve_atom_info(values[1], pdb)
                self.metrics.append(cs.DISTANCE_ATOMS_TAG.format(atom1, atom2, name))