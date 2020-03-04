import os
import shutil
import pele_platform.constants.constants as cs
import pele_platform.Utilities.Parameters.pele_env as pele
import pele_platform.Utilities.Helpers.simulation as ad
import pele_platform.Utilities.Helpers.constraints as cst
import pele_platform.Frag.helpers as hp
import pele_platform.Frag.fragments as fr

class FragRunner(pele.EnviroBuilder):

    #def __init__(self, core, input, gr_steps, steps, eq_steps):
    def __init__(self, args):
        self.software = "Frag"
        self.build_frag_variables(args)
        self.core = args.frag_core
        self.ligands = args.frag_ligands
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
            print(ligand.GetProp("_Name"))
            atom_core_idx = hp.search_core_fragment_linker(ligand, ligand_core)
            atom_name_core = ligand_core.GetAtomWithIdx(atom_core_idx).GetMonomerInfo().GetName()
            mol = Chem.MolFromPDBFile(self.core, removeHs=False)
            original = rd.SplitMolByPDBResidues(mol)[self.residue]
            for atom in original.GetAtomWithIdx(atom_core_idx).GetNeighbors(): print(atom.GetMonomerInfo().GetName())
            hydrogen_core_name = [atom.GetMonomerInfo().GetName() for atom in original.GetAtomWithIdx(atom_core_idx).GetNeighbors() if atom.GetAtomicNum() == 1][0]
            
            # Delete core for full ligand
            fragment = rd.DeleteSubstructs(ligand, ligand_core)
            Chem.MolToPDBFile(fragment, "int.pdb")
            
            new_mol = rc.EditableMol(fragment)
            for atom in reversed(fragment.GetAtoms()): 
                neighbours = atom.GetNeighbors()
                if len(neighbours) == 0: new_mol.RemoveAtom(atom.GetIdx())
            #Add missing hydrogen to full ligand and create pdb
            fragment = new_mol.GetMol()
            Chem.MolToPDBFile(fragment, "int2.pdb")
            old_atoms = [atom.GetIdx() for atom in fragment.GetAtoms()]
            fragment = Chem.AddHs(fragment, False, True)
            rp.EmbedMolecule(fragment)
            
            Chem.MolToPDBFile(fragment, "int3.pdb")
            added_hydrogen_idx = [atom.GetIdx() for atom in fragment.GetAtoms() if atom.GetIdx() not in old_atoms][0]
            atom_fragment_idx = fragment.GetAtomWithIdx(added_hydrogen_idx).GetNeighbors()[0].GetIdx()
            fragment_filename = fragment.GetProp("_Name")+".pdb"
            Chem.MolToPDBFile(fragment, fragment_filename)
            fragment = fr.Fragment(fragment_filename, atom_fragment_idx, added_hydrogen_idx, atom_name_core, hydrogen_core_name)
            self.fragments.append(fragment)
            lines.append(fragment.get_inputfile_line())
            fragment.santize_file()
        with open(self.input, "w") as fout:
            fout.write(("\n").join(lines))
        

    def set_test_variables(self):
        self.gr_steps = 1
        self.frag_steps = 1
        self.frag_eq_steps = 1
        self.temperature = 100000
        self.cpus = 4
	
