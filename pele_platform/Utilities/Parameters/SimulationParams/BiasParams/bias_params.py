
class BiasParams(object):


    def generate_bias_params(self, args):
        self.epsilon = args.epsilon if args.epsilon else self.simulation_params.get("epsilon", 0)
        self.bias_column = args.bias_column if args.bias_column else self.simulation_params.get("bias_column", 5)
