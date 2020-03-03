import os
import shutil
import pele_platform.constants.constants as cs
import pele_platform.Utilities.Parameters.pele_env as pele
import pele_platform.Utilities.Helpers.simulation as ad
import pele_platform.Utilities.Helpers.constraints as cst

class FragRunner(pele.EnviroBuilder):

    #def __init__(self, core, input, gr_steps, steps, eq_steps):
    def __init__(self, args):
        self.software = "Frag"
        self.build_frag_variables(args)
        self.core = args.frag_core
        self.input = args.frag_input
        self.gr_steps = args.growing_steps if args.growing_steps else self.simulation_params.get("growing_steps", 6)
        self.frag_steps = args.frag_steps if args.frag_steps else self.simulation_params.get("steps_in_gs", 3)
        self.frag_eq_steps = args.frag_eq_steps if args.frag_eq_steps else self.simulation_params.get("sampling_steps", 20)
        self.control_file = os.path.join(cs.DIR, "Templates/pele_template.conf")
        self.protocol = args.protocol if args.protocol else self.simulation_params.get("protocol", "")
        self.topology = None if self.pdb else os.path.join("output_pele", "topology.pdb")
        self.constraints = cst.retrieve_constraints(self.core, {}, {}, 5)
        self.box = cs.BOX.format(self.box_radius, self.box_center) if  self.box_radius else ""
        

    def prepare_control_file(self):
        tmp_control_file = os.path.basename(self.control_file)
        shutil.copy(self.control_file, tmp_control_file)
        adaptive = ad.SimulationBuilder("", tmp_control_file, self)
        # Fill to time because we have flags inside flags
        adaptive.fill_pele_template(self)
        adaptive.fill_pele_template(self)
        self.control_file = tmp_control_file

    def run(self):
        if self.protocol:
            command = "python -m frag_pele.main -cp {} -sef {} --sch_python {} --contrl {} -nc -d {} -dat {} -doc {} --license {} --cpus {} -{}".format(
                self.core, self.input, self.spython, self.control_file,
                self.pele_exec, self.pele_data, self.pele_documents, self.license,
                self.cpus, self.protocol)
        else:
            command = "python -m frag_pele.main -cp {} -sef {} --sch_python {} --contrl {} -nc -d {} -dat {} -doc {} --license {} --cpus {} --growing_steps {} --steps {} --pele_eq_steps {} --temperature  {}".format(
                self.core, self.input, self.spython, self.control_file,
                self.pele_exec, self.pele_data, self.pele_documents, self.license,
                self.cpus, self.gr_steps, self.frag_steps, self.frag_eq_steps, self.temperature)
        print(command)
        if not self.debug:
            os.system(command)

    def set_test_variables(self):
        self.gr_steps = 1
        self.frag_steps = 1
        self.frag_eq_steps = 1
        self.temperature = 100000
        self.cpus = 4
	
