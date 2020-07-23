import pele_platform.constants.constants as cs

class InOutParams(object):


    def generate_inout_params(self, args):
        #If User inputs exit condition or when doing an exit simulation
        if args.exit or args.in_out or args.in_out_soft:
            self.exit_value = args.exit_value if args.exit_value else self.simulation_params.get("exit_value", 0.9)
            self.exit_condition = args.exit_condition if args.exit_condition else self.simulation_params.get("exit_condition", ">")
            self.exit_trajnum = args.exit_trajnum if args.exit_trajnum else self.simulation_params.get("exit_trajnum", 4)
            self.bias_column = args.bias_column if args.bias_column else self.simulation_params.get("bias_column", 6)
            self.unbinding_block = cs.UNBINDING.format(self.bias_column, self.exit_value, self.exit_condition, self.exit_trajnum)
            self.equilibration = "false" #Lots of problems look into it
        else:
            self.unbinding_block = ""
