import pele_platform.Adaptive.metrics as mt


class FragMetrics(mt.MetricBuilder):


	def __init__(self, args):
            if args.atom_dist:
                self.metrics = self.distance_to_atom_json(self.system, args.atom_dist)
            else:
                self.metrics = ""
            if args.native:
                self.native = self.rsmd_to_json(args.native, self.chain)
            else:
                self.native = ""
