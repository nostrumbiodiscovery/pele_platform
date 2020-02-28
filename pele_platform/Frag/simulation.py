import os
import pele_platform.constants.constants as cs


class FragRunner():

    def __init__(self, core, input, cpus):
        self.core = core
        self.input = input
        self.control_file = os.path.join(cs.DIR, "Templates/frag_template.conf")
        self.spython = os.path.join(cs.SCHRODINGER, "utilities/python")
        if not os.path.exists(self.spython):
            self.spython = os.path.join(cs.SCHRODINGER, "run")
        self.pele_bin = cs.PELE_BIN
        self.data = os.path.join(cs.PELE, "Data")
        self.documents = os.path.join(cs.PELE, "Documents")
        self.license = cs.LICENSE
        self.cpus = cpus 
        self.test = True

    def run(self):
        command = "python -m frag_pele.main -cp {} -sef {} --sch_python {} --contrl {} -nc -d {} -dat {} -doc {} --license {} --cpus {}".format(
            self.core, self.input, self.spython, self.control_file, self.pele_bin, self.data, self.documents, self.license,
            self.cpus)
        if self.test:
            command += " --pele_eq_steps 1 --steps 1 --temp 1000000 --growing_steps 1 --restart"
        print(command)
        os.system(command)
	
