import pele_platform.Utilities.Helpers.helpers as hp
import pele_platform.constants.constants as cs

"""
Configuration of valid metrics for restrictions:
    - template: constant with json template for this metric.
    - description: string with a description about this metric.
    - num_elems: num of atoms required for this metric.
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
        """
        Parse interaction_restrictions (metrics and conditions) from configuration.

        Parameters
        ----------
        pdb : System pdb. Used to locate atoms.
        constraints_conf : Configuration parameters for interaction restrictions.
        """

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
        """
        Represent conditions for interaction restrictions in JSON format.

        Returns
        ----------
        json_string : str
            conditions formatted in JSON.
        """
        return cs.INTERACTION_RESTRICTIONS.format('",\n\t"'.join(self.conditions))

    def metrics_to_json(self):
        """
        Represent metrics for interaction restrictions in JSON format.

        Returns
        ----------
        json_string : str
            metrics formatted in JSON.
        """
        return "\n".join(self.metrics)

    def _add_metric(self, pdb, config, values, name):
        """
        Internal use only. Create a new metric.

        Parameters
        ----------
        pdb : System pdb. Used to locate atoms.
        config : Standard configuration for the metric.
        values : List of atoms.
        name : Name for the new metric.
        """
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
        """
        Internal use only. Create a new restriction condition.

        Parameters
        ----------
        values : Limits to apply in the condition (min and max).
        name : Name of the related metric.
        """
        if "min" in values:
            self.conditions.append(name + " > " + str(values["min"]))
        if "max" in values:
            self.conditions.append(name + " < " + str(values["max"]))

    def fill_template(self, template):
        """
        Joins self.conditions with AND and fills the PARAMETERS_CHANGE template string.

        Returns
        -------
        interaction_parameters_change : str
            Parameters change string for interaction restrictions, which will later be injected into pele_params.py.
        """
        joined_conditions = " and ".join(self.conditions) if len(self.conditions) > 1 else self.conditions[0]
        return template.format(joined_conditions)
