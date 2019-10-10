
class BiasParams(object):


    def __init__(self, args):
        self.epsilon = args.epsilon if args.epsilon else self.simulation_params.get("epsilon", None)
        self.bias_column = args.bias_column if args.bias_column else self.simulation_params.get("bias_column", None)
