import pele_platform.Utilities.Helpers.helpers as hp
import pele_platform.constants.constants as cs

"""
Configuration of valid metrics for restrictions
"""
RESTRICTIONS_CONFIG = {
    "distance": {
        "template": cs.DISTANCE_ATOMS_TAG,
        "description": "distance",
        "num_elems": 2,
    },
    "angle": {"template": cs.ANGLE_ATOMS_TAG, "description": "angle", "num_elems": 3},
}


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
        for i, constraint_conf in enumerate(constraints_conf):
            if constraint_conf:
                restriction = set(RESTRICTIONS_CONFIG.keys()).intersection(
                    constraint_conf.keys()
                )
            if len(restriction) == 1:
                id = restriction.pop()
                name = id + str(i)
                self._add_metric(
                    pdb, RESTRICTIONS_CONFIG[id], constraint_conf["atoms"], name
                )
                self._create_conditions(constraint_conf[id], name)
            else:
                raise SyntaxError(
                    f"Must define only one valid metric in each restriction. The valid "
                    f"metrics are: {list(RESTRICTIONS_CONFIG)}. You defined {len(restriction)} "
                    f"in {constraint_conf}."
                )

    def conditions_to_json(self):
        return cs.INTERACTION_RESTRICTIONS.format('",\n\t"'.join(self.conditions))

    def metrics_to_json(self):
        return "\n".join(self.metrics)

    def _add_metric(self, pdb, config, values, name):
        num_atoms = config["num_elems"]
        if len(values) != num_atoms:
            raise SyntaxError(
                f"Must specify a list of {num_atoms} atoms in {config['description']} restriction."
            )
        atoms = []
        for atom in values:
            atoms.append(hp.retrieve_atom_info(atom, pdb))
        self.metrics.append(config["template"].format(name, *atoms))

    def _create_conditions(self, values, name):
        if "min" in values:
            self.conditions.append(name + " > " + str(values["min"]))
        if "max" in values:
            self.conditions.append(name + " < " + str(values["max"]))
