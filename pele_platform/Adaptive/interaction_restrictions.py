import pele_platform.Utilities.Helpers.helpers as hp
import pele_platform.constants.constants as cs
import os


class InteractionRestrictionsBuilder():

    def __init__(self):
        self.metrics = []
        self.conditions = []

    def parse_interaction_restrictions (self, pdb, constraints_conf):
        for i in range (0, len(constraints_conf)):
            actual = constraints_conf[i]
            name = ""
            if 'distance' in actual:
                name = "distance" + str(i)
                self._add_metric (pdb, "atom_dist", actual['atoms'], name)
                self._create_conditions (actual['distance'], name)
            if 'angle' in actual:
                name = "angle" + str(i)
                self._add_metric (pdb, "atom_angle", actual['atoms'], name)
                self._create_conditions (actual['angle'], name)

    def conditions_to_json (self):
        return cs.INTERACTION_RESTRICTIONS.format('",\n\t"'.join(self.conditions))

    def metrics_to_json (self):
        return "\n".join(self.metrics)

    def _add_metric (self, pdb, type, values, name):
        if (type == "atom_dist"):
                atom1 = hp.retrieve_atom_info(values[0], pdb)
                atom2 = hp.retrieve_atom_info(values[1], pdb)
                self.metrics.append(cs.DISTANCE_ATOMS_TAG.format(atom1, atom2, name))
        if (type == "atom_angle"):
                atom1 = hp.retrieve_atom_info(values[0], pdb)
                atom2 = hp.retrieve_atom_info(values[1], pdb)
                atom3 = hp.retrieve_atom_info(values[2], pdb)
                self.metrics.append(cs.ANGLE_ATOMS_TAG.format(atom1, atom2, atom3, name))

    def _create_conditions (self, values, name):
        if 'min' in values:
            self.conditions.append(name + " > " + str (values['min']))
        if 'max' in values:
            self.conditions.append(name + " < " + str (values['max']))
