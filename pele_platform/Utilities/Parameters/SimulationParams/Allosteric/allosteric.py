
class AllostericParams(object):


    def generate_allosteric_params(self, args):
        self.allosteric = args.allosteric
        self.n_components = args.n_components if args.n_components else  self.simulation_params.get("n_components", 10)
