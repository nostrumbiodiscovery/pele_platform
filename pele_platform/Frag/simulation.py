import os
import tempfile
import shutil
import pele_platform.Utilities.Helpers.simulation as ad
import pele_platform.Frag.helpers as hp
import pele_platform.Frag.checker as ch
import pele_platform.Frag.parameters.main as mn
import pele_platform.Errors.custom_errors as ce
import frag_pele.main as frag


class FragRunner(mn.FragParameters):

    def __init__(self, args):
        args.system = args.frag_core
        mn.FragParameters.__init__(self, args)

    def run_simulation(self):
        self._set_test_variables()
        self._prepare_control_file()
        self._launch()

    def _launch(self):
        if self.ligands:  # Full ligands as sdf
            fragment_files = self._prepare_input_file(logger=self.logger)
        else:
            fragment_files = None
        self._run()
        
        if self.cleanup and fragment_files:
            self._clean_up(fragment_files)

    def _prepare_control_file(self):
        # Create tmp folder with frag control_file 
        tmp_dir = tempfile.mkdtemp()
        tmp_control_file = os.path.join(tmp_dir, os.path.basename(self.control_file))
        shutil.copy(self.control_file, tmp_control_file)
        adaptive = ad.SimulationBuilder("", tmp_control_file, self)
        # Fill to time because we have flags inside flags
        adaptive.fill_pele_template(self, self.water_object)
        adaptive.fill_pele_template(self, self.water_object)
        self.control_file = tmp_control_file
        return self.control_file

    def _set_test_variables(self):
        if self.test:
            self.gr_steps = 1
            self.frag_steps = 2
            self.frag_eq_steps = 1
            self.temperature = 100000
            self.anm_freq = 0
            self.minimizatoon = 0
            self.sidechain_freq = 0
            self.water_freq = 0
            self.cpus = 4
        else:
            pass

    def _run(self):
        if self.frag_run:
            try:
                frag.main(self.core_process, self.input, self.gr_steps, self.criteria, self.plop_path, self.spython,
                          self.pele_exec, self.control_file, self.license, self.output_folder,
                          self.report_name, "trajectory", self.cluster_folder, self.cpus, self.distcont, self.threshold,
                          self.epsilon, self.condition, self.metricweights,
                          self.nclusters, self.frag_eq_steps, self.frag_restart, self.min_overlap, self.max_overlap,
                          self.chain_core, self.frag_chain, self.frag_steps, self.temperature, self.seed, self.gridres,
                          self.banned, self.limit, self.mae,
                          self.rename, self.threshold_clash, self.steering, self.translation_high, self.rotation_high,
                          self.translation_low, self.rotation_low, self.explorative, self.frag_radius,
                          self.sampling_control, self.pele_data, self.pele_documents,
                          self.only_prepare, self.only_grow, self.no_check, self.debug, srun=self.usesrun)
            except Exception:
                print("Skipped - FragPELE will not run.")

    def _prepare_input_file(self, logger=None):
        from rdkit import Chem
        import rdkit.Chem.rdmolops as rd

        self.input = "input.conf"

        # Check input file
        limit_atoms_ligand = 100
        ch.check_limit_number_atoms(self.ligands, limit_atoms_ligand)
        
        # Get core of the ligand
        mol = Chem.MolFromPDBFile(self.core)
        try:
            ligand_core = rd.SplitMolByPDBResidues(mol)[self.residue]
        except KeyError:
            raise ce.MissResidueFlag("Missing residue flag to specify the ligand core residue name. i.e resname: 'LIG'")
            

        # Get sdf full grown ligands
        ligands_grown = Chem.SDMolSupplier(self.ligands, removeHs=False)
        fragment_files = []
        
        # For each full grown ligand create neutral fragment
        with open(self.input, "w") as fout:
            pass
        for ligand in ligands_grown:
            Chem.AssignAtomChiralTagsFromStructure(ligand)
            try:
                line, fragment = self._create_fragment_from_ligand(ligand, ligand_core, logger=logger)
            except Exception as e:
                try:
                    line, fragment = self._create_fragment_from_ligand(ligand, ligand_core, substructure=False)
                except Exception as e:
                    try:
                        # Try to fix symmetry
                        line, fragment = self._create_fragment_from_ligand(ligand, ligand_core, symmetry=True)
                    except Exception as e:
                            # Try with second substructure search
                        line, fragment = self._create_fragment_from_ligand(ligand, ligand_core, result=1, substructure=False)

            logger.info(f"Ligand {fragment.file} preprocessed")
            fragment_files.append(fragment.file)
            with open(self.input, "a") as fout:
                fout.write(line + "\n")

        return fragment_files

    def _create_fragment_from_ligand(self, ligand, ligand_core, result=0, substructure=True, symmetry=False):
        import rdkit.Chem.AllChem as rp
        from rdkit import Chem

        fragment, old_atoms, hydrogen_core, atom_core, atom_frag, mapping, correct = hp._build_fragment_from_complex(
            self.core, self.residue, ligand, ligand_core, result, substructure, symmetry)

        # temporary override to fix segmentation faults
        filename = "temp.pdb"
        Chem.MolToPDBFile(fragment, filename)
        fragment = Chem.MolFromPDBFile(filename, removeHs=False)
        os.remove(filename)
        
        rp.EmbedMolecule(fragment)
        fragment = hp._retrieve_fragment(
        fragment, old_atoms, atom_core, hydrogen_core, atom_frag, mapping)
        line = fragment.get_inputfile_line()
        fragment.sanitize_file()
        
        if not correct:
            print("Ligand incorrect")
        return line, fragment
        

    def _clean_up(self, fragment_files):

        for file in fragment_files:
            if os.path.isfile(file):
                os.remove(file)
