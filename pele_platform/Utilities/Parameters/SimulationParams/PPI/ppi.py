
class PPIParams(object):


    def __init__(self, args):
        self.ppi = args.ppi
        self.n_components = args.n_components if args.n_components else  self.simulation_params.get("n_components", 10)
