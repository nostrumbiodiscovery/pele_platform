import os
import tempfile
import shutil
import pele_platform.Utilities.Helpers.simulation as ad
import pele_platform.Frag.helpers as hp
import pele_platform.Frag.fragments as fr
import pele_platform.Frag.checker as ch
import pele_platform.Frag.generative_model as gm
import pele_platform.Frag.parameters.main as mn
import pele_platform.Errors.custom_errors as ce
import frag_pele.main as frag

class FragRunner(mn.FragParameters):

    #def __init__(self, core, input, gr_steps, steps, eq_steps):
    def __init__(self, args):
        mn.FragParameters.__init__(self, args)

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
         if self.frag_run:
             frag.main(self.core_process, self.input, self.gr_steps, self.criteria, self.plop_path, self.spython, self.pele_exec, self.control_file, self.license, self.output_folder,
             self.report_name, "trajectory", self.cluster_folder, self.cpus, self.distcont, self.threshold, self.epsilon, self.condition, self.metricweights,
             self.nclusters, self.frag_eq_steps, self.frag_restart, self.min_overlap, self.max_overlap,
             self.chain_core, self.frag_chain, self.frag_steps, self.temperature, self.seed, self.gridres, self.banned, self.limit, self.mae,
             self.rename, self.threshold_clash, self.steering, self.translation_high, self.rotation_high,
             self.translation_low, self.rotation_low, self.explorative, self.frag_radius, self.sampling_control, self.pele_data, self.pele_documents,
             self.only_prepare, self.only_grow, self.no_check, self.debug)

    def prepare_input_file(self):
        from rdkit import Chem
        import rdkit.Chem.rdmolops as rd
        import rdkit.Chem.rdchem as rc
        import rdkit.Chem.AllChem as rp

        self.input = "input.conf"

        #Check input file
        limit_atoms_ligand = 100
        ch.check_limit_number_atoms(self.ligands, limit_atoms_ligand)
        
        #Get core of the ligand
        mol = Chem.MolFromPDBFile(self.core)
        try:
            ligand_core = rd.SplitMolByPDBResidues(mol)[self.residue]
        except KeyError:
            raise ce.MissResidueFlag("Missing residue flag to specify the ligand core residue name. i.e resname: 'LIG'")
            
        #Get sdf full grown ligands
        ligands_grown = Chem.SDMolSupplier(self.ligands, removeHs=False)

        #For each full grown ligand create neutral fragment
        self.fragments = []
        lines = []
        for ligand in ligands_grown:
            try:
                line, fragment = self._create_fragment_from_ligand(ligand, ligand_core)
            except Exception as e:
                try:
                    # Try to fix simmetry
                    line, fragment = self._create_fragment_from_ligand(ligand, ligand_core, simmetry=True)
                except Exception as e:
                    try:
                        # Try with second substructure search
                        line, fragment = self._create_fragment_from_ligand(ligand, ligand_core, result=1, substructure=False)
                    except Exception as e:
                        #Skip the ligand
                        print("Ligand Skipped")
                        print(e)
                        continue
            lines.append(line)
            self.fragments.append(fragment)
            print(f"Ligand {fragment.file} preprocessed")
        with open(self.input, "w") as fout:
            fout.write(("\n").join(lines))


    def _create_fragment_from_ligand(self, ligand, ligand_core, result=0, substructure=True, simmetry=False):
        from rdkit import Chem
        import rdkit.Chem.rdmolops as rd
        import rdkit.Chem.rdchem as rc
        import rdkit.Chem.AllChem as rp
        fragment, old_atoms, hydrogen_core, atom_core = hp._build_fragment_from_complex(
            self.core, self.residue, ligand, ligand_core, result, substructure, simmetry)
        rp.EmbedMolecule(fragment)
        fragment = hp._retrieve_fragment(
            fragment, old_atoms, atom_core, hydrogen_core)
        line = fragment.get_inputfile_line()
        fragment.sanitize_file()
        return line, fragment
        

    def set_test_variables(self):
        self.gr_steps = 1
        self.frag_steps = 2
        self.frag_eq_steps = 1
        self.temperature = 100000
        self.cpus = 4
	
