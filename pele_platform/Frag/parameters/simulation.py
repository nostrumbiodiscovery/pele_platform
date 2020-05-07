import os
import pele_platform.constants.constants as cs

class FragSimulationParameters():

    def __init__(self, args):
        self.gr_steps = args.growing_steps if args.growing_steps else self.simulation_params.get("growing_steps", 6)
        self.frag_steps = args.frag_steps if args.frag_steps else self.simulation_params.get("steps_in_gs", 3)
        self.frag_eq_steps = args.frag_eq_steps if args.frag_eq_steps else self.simulation_params.get("sampling_steps", 20)
        self.control_file = os.path.join(cs.DIR, "Templates/pele_template.conf")
        self.protocol = args.protocol if args.protocol else self.simulation_params.get("protocol", "")
        self.topology = None if self.pdb else os.path.join("output_pele", "topology.pdb")
        self.frag_restart = "-rst" if args.frag_restart else ""
