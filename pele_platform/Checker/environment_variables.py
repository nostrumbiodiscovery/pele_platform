import os

class EnvVariable():

    def __init__(self, name, variable, default, flag, env_var):
        self.name = name
        self.variable = variable
        self.default = default
        self.flag = flag
        self.env_var = env_var

    def is_valid(self):
        if self.variable:
            if os.path.exists(self.variable):
                return True
        elif self.default:
           if os.path.exists(self.default):
               return True
        raise ValueError("{} not found. \n\
1) If you have the standard installation export the environment variable by doing:\n\
\t{}.\n\
2) If the previous does not work define the location of {} with the next flag in the input.yaml:\n\
\t{}".format(self.name, self.env_var, self.name, self.flag))

