import os
import pele_platform.constants.constants as cs
import pele_platform.Utilities.Parameters.pele_env as pele


class FragRunner(pele.EnviroBuilder):

    #def __init__(self, core, input, gr_steps, steps, eq_steps):
    def __init__(self, args):
        self.software = "Frag"
        self.build_frag_variables(args)
        self.core = args.frag_core
        self.input = args.frag_input
        self.gr_steps = args.growing_steps if args.growing_steps else self.simulation_params.get("growing_steps", 6)
        self.steps = args.steps if args.steps else self.simulation_params.get("steps", 3)
        self.eq_steps = args.eq_steps if args.eq_steps else self.simulation_params.get("eq_steps", 20)
        self.control_file = os.path.join(cs.DIR, "Templates/frag_template.conf")

    def run(self):
        command = "python -m frag_pele.main -cp {} -sef {} --sch_python {} --contrl {} -nc -d {} -dat {} -doc {} --license {} --cpus {} --growing_steps {} --steps {} --pele_eq_steps {} --temperature  {}".format(
            self.core, self.input, self.spython, self.control_file, self.pele_exec, self.pele_data, self.pele_documents, self.license,
            self.cpus, self.gr_steps, self.steps, self.eq_steps, self.temperature)
        print(command)
        os.system(command)

    def set_test_variables(self):
        self.gr_steps = 1
        self.steps = 1
        self.eq_steps = 1
        self.temperature = 10000
        self.cpus = 2
	
