import pele_platform.Utilities.Helpers.helpers as hp
import pele_platform.constants.constants as cs

RESTRICTIONS_CONFIG = {
    "distance": {
        "template": cs.DISTANCE_ATOMS_TAG,
        "description": "distance",
        "num_elems": 2,
    },
    "angle": {"template": cs.ANGLE_ATOMS_TAG, "description": "angle", "num_elems": 3},
}
"""
Configuration of valid metrics for restrictions
"""


class InteractionRestrictionsBuilder:
    """
    Base class to generate interaction restrictions
    1) Parse interaction restrictions block (metrics and conditions)
    2) Generate json format for pele_template.conf
    """

    def __init__(self):
        self.metrics = []
        self.conditions = []

    def parse_interaction_restrictions(self, pdb, constraints_conf):
        for i in range(0, len(constraints_conf)):
            actual = constraints_conf[i]
            restriction = set(RESTRICTIONS_CONFIG.keys()).intersection(actual.keys())
            if len(restriction) == 1:
                id = restriction.pop()
                name = id + str(i)
                self._add_metric(pdb, RESTRICTIONS_CONFIG[id], actual["atoms"], name)
                self._create_conditions(actual[id], name)

    def conditions_to_json(self):
        return cs.INTERACTION_RESTRICTIONS.format('",\n\t"'.join(self.conditions))

    def metrics_to_json(self):
        return "\n".join(self.metrics)

    def _add_metric(self, pdb, config, values, name):
        num_atoms = config["num_elems"]
        if len(values) != num_atoms:
            raise SyntaxError(
                "Must specify a list of "
                + str(num_atoms)
                + " atoms in "
                + config["description"]
                + " restriction."
            )
        atoms = []
        for i in range(num_atoms):
            atoms.append(hp.retrieve_atom_info(values[i], pdb))
        self.metrics.append(config["template"].format(name, *atoms))

    def _create_conditions(self, values, name):
        if "min" in values:
            self.conditions.append(name + " > " + str(values["min"]))
        if "max" in values:
            self.conditions.append(name + " < " + str(values["max"]))
