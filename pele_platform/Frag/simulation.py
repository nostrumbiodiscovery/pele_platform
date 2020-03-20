import os
import tempfile
import shutil
import pele_platform.constants.constants as cs
import pele_platform.Utilities.Parameters.pele_env as pele
import pele_platform.Utilities.Helpers.simulation as ad
import pele_platform.Utilities.Helpers.constraints as cst
import pele_platform.Frag.helpers as hp
import pele_platform.Frag.fragments as fr
import pele_platform.Frag.generative_model as gm

class FragRunner(pele.EnviroBuilder):

    #def __init__(self, core, input, gr_steps, steps, eq_steps):
    def __init__(self, args):
        self.software = "Frag"
        self.build_frag_variables(args)
        self.core = args.frag_core
        self.ligands = args.frag_ligands
        self.ai = args.frag_ai
        self.ai_iterations = args.frag_ai_iterations
        self.core_format = args.frag_core.rsplit(".")[-1]
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
        # Create tmp folder with frag control_file 
        tmp_dir = tempfile.mkdtemp()
        tmp_control_file = os.path.join(tmp_dir, os.path.basename(self.control_file))
        shutil.copy(self.control_file, tmp_control_file)
        adaptive = ad.SimulationBuilder("", tmp_control_file, self)
        # Fill to time because we have flags inside flags
        adaptive.fill_pele_template(self)
        adaptive.fill_pele_template(self)
        self.control_file = tmp_control_file
        return self.control_file

    def run(self):
        # If protocol let frag handle all flags
        if self.protocol:
            command = "python -m frag_pele.main -cp {} -sef {} --sch_python {} --contrl {} -d {} -dat {} -doc {} --license {} --cpus {} -{}".format(
                self.core, self.input, self.spython, self.control_file,
                self.pele_exec, self.pele_data, self.pele_documents, self.license,
                self.cpus, self.protocol)
        else:
            # Pass all possible flags
            command = "python -m frag_pele.main -cp {} -sef {} --sch_python {} --contrl {} -d {} -dat {} -doc {} --license {} --cpus {} --growing_steps {} --steps {} --pele_eq_steps {} --temperature  {}".format(
                self.core, self.input, self.spython, self.control_file,
                self.pele_exec, self.pele_data, self.pele_documents, self.license,
                self.cpus, self.gr_steps, self.frag_steps, self.frag_eq_steps, self.temperature)
        print(command)
        if not self.debug:
            os.system(command)

    def prepare_input_file(self):
        from rdkit import Chem
        import rdkit.Chem.rdmolops as rd
        import rdkit.Chem.rdchem as rc
        import rdkit.Chem.AllChem as rp

        self.input = "input.conf"
        
        #Get core of the ligand
        mol = Chem.MolFromPDBFile(self.core)
        ligand_core = rd.SplitMolByPDBResidues(mol)[self.residue]
            
        #Get sdf full grown ligands
        ligands_grown = Chem.SDMolSupplier(self.ligands, removeHs=False)

        #For each full grown ligand create neutral fragment
        self.fragments = []
        lines = []
        for ligand in ligands_grown:
            try:
                line, fragment = self._create_fragment_from_ligand(ligand, ligand_core)
            except (ValueError, AttributeError):
                try:
                    line, fragment = self._create_fragment_from_ligand(ligand, ligand_core, result=1, substructure=False)
                except (IndexError):
                    print("LIGAND SKIPPED")
                    continue
            
            lines.append(line)
            self.fragments.append(fragment)
        with open(self.input, "w") as fout:
            fout.write(("\n").join(lines))


    def _create_fragment_from_ligand(self, ligand, ligand_core, result=0, substructure=True):
        from rdkit import Chem
        import rdkit.Chem.rdmolops as rd
        import rdkit.Chem.rdchem as rc
        import rdkit.Chem.AllChem as rp
        fragment, old_atoms, hydrogen_core, atom_core = hp._build_fragment_from_complex(
            self.core, self.residue, ligand, ligand_core, result, substructure)
        rp.EmbedMolecule(fragment)
        fragment = hp._retrieve_fragment(
            fragment, old_atoms, atom_core, hydrogen_core)
        line = fragment.get_inputfile_line()
        fragment.sanitize_file()
        return line, fragment
        

    def grow_ai(self):
        with open(self.core, "r") as fin:
            ligand_lines = [line for line in fin if line[17:20] == self.residue or line.startswith("CONECT")]
            ligand_lines.append("END")
        self.ai = "core.pdb"
        with open(self.ai, "w") as fout:
            fout.write("".join(ligand_lines))
            
           
        generative_model_rnn = gm.GenerativeModel(self.ai, self.residue, self.ai_iterations)
        ligands = generative_model_rnn.run()
        for ligand in ligands:
            self.ligands = ligand
            self.prepare_input_file()
            self.run()

    def set_test_variables(self):
        self.gr_steps = 1
        self.frag_steps = 1
        self.frag_eq_steps = 1
        self.temperature = 100000
        self.cpus = 4
	
