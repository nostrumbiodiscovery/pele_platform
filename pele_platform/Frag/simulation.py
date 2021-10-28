import os
import tempfile
import shutil
import re
import glob


import pele_platform.Utilities.Helpers.simulation as ad
import pele_platform.Frag.helpers as hp
import pele_platform.Frag.checker as ch
import pele_platform.Errors.custom_errors as ce
import pele_platform.Frag.libraries as lb
import pele_platform.Frag.analysis as ana
import pele_platform.constants.constants as cs
import pele_platform.Frag.filtering as fl
import frag_pele.main as frag
from pele_platform.analysis import Analysis
from pele_platform.Utilities.Helpers import water



class FragRunner(object):
    def __init__(self, args):
        from pele_platform.Utilities.Parameters import ParametersBuilder

        builder = ParametersBuilder()
        parameters = builder.build_frag_variables(args)
        self._parameters = parameters

    @property
    def parameters(self):
        return self._parameters

    def run_simulation(self):
        self._set_test_variables()
        self._prepare_control_file()
        if not self.parameters.only_filtering:
            self._launch()
            if not self.parameters.debug:
                self._analysis()
        if self.parameters.database:
            self._filtering()
        return self.parameters

    def _launch(self):
        params = self.parameters
        with tempfile.TemporaryDirectory() as tmpdirname:
            if params.ligands:  # Full ligands as sdf
                fragment_files = self._prepare_input_file(logger=params.logger)
            elif params.frag_library:
                params.input = lb.main(params.frag_core_atom,
                                       params.frag_library,
                                       params.logger,
                                       params.fragment_atom,
                                       tmpdirname)
                if params.frag_restart_libraries:
                    self._frag_restart()

            else:
                fragment_files = None
                assert params.input is not None, (
                    "You need to provide input.conf file or SD file with fully grown "
                    "ligands. "
                )

            if not params.only_analysis:
                self._prepare_parameters()
                self._run()

            if params.cleanup and fragment_files:
                self._clean_up(fragment_files)

    def _prepare_control_file(self):
        # Create tmp folder with frag control_file
        params = self.parameters

        tmp_dir = tempfile.mkdtemp()
        tmp_control_file = os.path.join(tmp_dir, os.path.basename(params.control_file))
        shutil.copy(params.control_file, tmp_control_file)
        adaptive = ad.SimulationBuilder("", tmp_control_file, params)
        # Fill to time because we have flags inside flags
        adaptive.fill_pele_template(params, params.water_object)
        adaptive.fill_pele_template(params, params.water_object)
        self.parameters.control_file = tmp_control_file

        return self.parameters.control_file

    def _prepare_parameters(self):
        self.parameters.spython = cs.SCHRODINGER  # Commented to use Frag 2.2.1 instead of Frag 3.0.0

    def _set_test_variables(self):
        if self.parameters.test:
            self.parameters.gr_steps = 1
            self.parameters.frag_steps = 2
            self.parameters.frag_eq_steps = 1
            self.parameters.temperature = 100000
            self.parameters.anm_freq = 0
            self.parameters.minimizatoon = 0
            self.parameters.sidechain_freq = 0
            self.parameters.water_freq = 0
            self.parameters.cpus = 4
        else:
            pass

    def _run(self):
        params = self.parameters
        if not params.frag_restart_libraries:
            self._extract_working_directory()
        try:
            frag.main(
                params.core_process,
                params.input,
                params.gr_steps,
                params.criteria,
                params.plop_path,
                params.spython,
                params.pele_exec,
                params.control_file,
                params.license,
                params.output_folder,
                params.report_name,
                "trajectory",
                params.cluster_folder,
                params.cpus,
                params.distcont,
                params.threshold,
                params.epsilon,
                params.condition,
                params.metricweights,
                params.nclusters,
                params.frag_eq_steps,
                params.frag_restart,
                params.min_overlap,
                params.max_overlap,
                params.chain_core,
                params.frag_chain,
                params.frag_steps,
                params.temperature,
                params.seed,
                params.gridres,
                params.banned,
                params.limit,
                params.mae,
                params.rename,
                params.threshold_clash,
                params.steering,
                params.translation_high,
                params.rotation_high,
                params.translation_low,
                params.rotation_low,
                params.explorative,
                params.frag_radius,
                params.sampling_control,
                params.pele_data,
                params.pele_documents,
                params.only_prepare,
                params.only_grow,
                params.no_check,
                params.debug,
                srun=params.usesrun,
            )
        except Exception as e:
            print("Skipped - FragPELE will not run.")
            print(e)

    def _prepare_input_file(self, logger=None):
        from rdkit import Chem
        import rdkit.Chem.rdmolops as rd

        self.parameters.input = "input.conf"

        params = self.parameters

        # Check input file
        limit_atoms_ligand = 100
        ch.check_limit_number_atoms(params.ligands, limit_atoms_ligand)

        # Get core of the ligand
        mol = Chem.MolFromPDBFile(params.core)
        try:
            ligand_core = rd.SplitMolByPDBResidues(mol)[params.residue]
        except KeyError:
            raise ce.MissResidueFlag(
                "Missing residue flag to specify "
                + "the ligand core residue name. "
                + "i.e resname: 'LIG'"
            )

        # Get sdf full grown ligands
        ligands_grown = Chem.SDMolSupplier(params.ligands, removeHs=False)
        fragment_files = []

        # For each full grown ligand create neutral fragment
        with open(params.input, "w") as fout:
            pass
        for ligand in ligands_grown:
            Chem.AssignAtomChiralTagsFromStructure(ligand)
            try:
                line, fragment = self._create_fragment_from_ligand(ligand, ligand_core)
            except Exception as e:
                try:
                    line, fragment = self._create_fragment_from_ligand(
                        ligand, ligand_core, substructure=False
                    )
                except Exception as e:
                    try:
                        # Try to fix symmetry
                        line, fragment = self._create_fragment_from_ligand(
                            ligand, ligand_core, symmetry=True
                        )
                    except Exception as e:
                        # Try with second substructure search
                        line, fragment = self._create_fragment_from_ligand(
                            ligand, ligand_core, result=1, substructure=False
                        )

            logger.info(f"Ligand {fragment.file} preprocessed")
            fragment_files.append(fragment.file)
            with open(params.input, "a") as fout:
                fout.write(line + "\n")

        return fragment_files

    def _create_fragment_from_ligand(
        self, ligand, ligand_core, result=0, substructure=True, symmetry=False
    ):
        import rdkit.Chem.AllChem as rp
        from rdkit import Chem

        params = self.parameters
        (
            fragment,
            old_atoms,
            hydrogen_core,
            atom_core,
            atom_frag,
            mapping,
            correct,
        ) = hp._build_fragment_from_complex(
            params.core,
            params.residue,
            ligand,
            ligand_core,
            result,
            substructure,
            symmetry,
        )

        # temporary override to fix segmentation faults
        filename = "temp.pdb"
        Chem.MolToPDBFile(fragment, filename)
        fragment = Chem.MolFromPDBFile(filename, removeHs=False)
        os.remove(filename)

        rp.EmbedMolecule(fragment)
        fragment = hp._retrieve_fragment(
            fragment, old_atoms, atom_core, hydrogen_core, atom_frag, mapping
        )
        line = fragment.get_inputfile_line()
        fragment.sanitize_file()

        if not correct:
            print("Ligand incorrect")
        return line, fragment

    def _analysis(self):
        # TODO create a list of the libraries defined in the current input.yaml
        # Retrieve water indices to cluster, if running analysis only
        if self.parameters.only_analysis:
            self._extract_working_directory()
            for path in self.parameters.working_dir:
                if not os.path.exists(path):
                    if not os.path.exists(self.parameters.folder):
                        raise ce.MissingWorkingDir("Missing working_folder parameter. Please set the "
                                                   "path of the trajectories using the flag working_"
                                                   "folder in your input.yaml")
            self._prepare_parameters()
            self.parameters.water_ids_to_track = water.water_ids_from_conf(self.parameters.control_file)

        if self.parameters.args.analysis_to_point:
            self.parameters.analysis_to_point = self.parameters.args.analysis_to_point
            ana.main(
                path=self.parameters.folder,
                atom_coords=self.parameters.analysis_to_point,
                pattern=os.path.basename(self.parameters.system),
            )
        for path in self.parameters.working_dir:
            simulation_output = os.path.join(path, 'sampling_result')
            analysis_folder = os.path.join(path, "results")

            analysis = Analysis(
                resname="GRW",
                chain=self.parameters.chain,
                simulation_output=simulation_output,
                be_column=self.parameters.be_column,
                limit_column=self.parameters.limit_column,
                traj=self.parameters.traj_name,
                report=self.parameters.report_name,
                skip_initial_structures=not self.parameters.test,
                kde=self.parameters.kde,
                kde_structs=self.parameters.kde_structs,
                topology=self.parameters.topology,
                cpus=self.parameters.cpus,
                water_ids_to_track=self.parameters.water_ids_to_track,
            )
            analysis.generate(
                analysis_folder,
                clustering_type=self.parameters.clustering_method.lower(),
                bandwidth=self.parameters.bandwidth,
                analysis_nclust=self.parameters.analysis_nclust,
                max_top_clusters=self.parameters.max_top_clusters,
                top_clusters_criterion=self.parameters.top_clusters_criterion,
                min_population=self.parameters.min_population,
                representatives_criterion=self.parameters.cluster_representatives_criterion,
            )

    def _filtering(self):
        """
        Runs Feed Forward Frag algorithm.
        The F3 algorithm can run in two ways:
            1. Performing a fragment growing simulation before running F3 algorithm.
            2. Running F3 algorithm directly on a protein-ligand system.
        Here we extract the initial system and the ligand to then start the filtering
        process where we obtain a list with the filtering results.
        For each molecule retrieved from the filtering process, the algorithm performs
        a growing simulation.
        If the molecule has more than one R-group bound to the core (our intitial seed
        compound), several fragment growing simulations are performed in order to grow
        all the R-groups. To do so, it is needed to change the directory of the PDB with
        the system as more R-groups are being bound to the core.
        """
        from rdkit import Chem
        top_results = []
        binding_energies = {}
        initial_system_min_energy = self.parameters.core
        initial_ligand_min_energy = self.parameters.f3_ligand
        self.parameters.input = "input.conf"
        initial_ligand_min_energy_mol = Chem.rdmolfiles.MolFromPDBFile(initial_ligand_min_energy, removeHs= True)
        Chem.rdmolops.AssignStereochemistryFrom3D(initial_ligand_min_energy_mol)
        top_molecules = fl.main(initial_ligand_min_energy_mol, initial_ligand_min_energy,
                               self.parameters.database, self.parameters.filters)
        print("You will perform a total of %s simulations" % len(top_molecules))
        print("System min energy:", initial_system_min_energy)
        print("ligand min energy:", initial_ligand_min_energy)
        for molecule in top_molecules:
            for fragment in range(len(top_molecules[molecule])):
                if fragment > 0:
                    self._extract_working_directory()
                    for dir in self.parameters.working_dir:
                        for result in glob.glob(dir + "/top_result/*"):
                            top_results.append(result)
                    for top_result in top_results:
                        binding_energies[top_result] = os.path.splitext(top_result)[0].split('y')[-1]
                    system_min_energy = min(binding_energies, key=binding_energies.get)
                    ligand_min_energy = system_min_energy.split("/")[:-2][-1] + '/pregrow/growing_result_p.pdb'
                    ligand_min_energy_mol = Chem.rdmolfiles.MolFromPDBFile(ligand_min_energy, removeHs=True)
                    Chem.rdmolops.AssignStereochemistryFrom3D(ligand_min_energy_mol)
                    self.parameters.core = system_min_energy
                    self.parameters.core_process = system_min_energy
                else:
                    self.parameters.core = initial_system_min_energy
                    self.parameters.core_process = initial_system_min_energy

                self.parameters.frag_library = os.path.join("./frag_library", molecule + "_" + str(fragment))
                self.parameters.frag_core_atom = top_molecules[molecule][fragment][0]
                self.parameters.fragment_atom = top_molecules[molecule][fragment][1]
                self._prepare_input_file_filtering()
                self._launch()

    def _clean_up(self, fragment_files):
        for file in fragment_files:
            if os.path.isfile(file):
                os.remove(file)

    def _prepare_input_file_filtering(self):
        with open("input.conf", "w+") as conf_file:
            conf_file.write(self.parameters.core + " " + self.parameters.frag_core_atom + " " + self.parameters.fragment_atom + "\n")

    def _extract_working_directory(self):
        """
        Extracts current working directory with FragPELE simulation's output.
        """
        params = self.parameters
        params.working_dir = []
        if os.path.isfile(params.core):
            complex_name = os.path.basename(params.core).split(".pdb")[0]  # And if it is a path, get only the name
        else:
            complex_name = params.core.split(".pdb")[0]
        pdb_basename = complex_name
        current_path = os.path.abspath(".")
        with open(params.input, "r") as input_file:
            for line in input_file.readlines():
                ID = os.path.basename(line).replace(".pdb", "")
                sentence = re.sub(r"\s+", "", ID, flags=re.UNICODE)
                if self.parameters.folder:
                    params.working_dir.append(os.path.join(self.parameters.folder, "{}_{}".format(pdb_basename, sentence)).strip('\n'))
                else:
                    params.working_dir.append(os.path.join(current_path, "{}_{}".format(pdb_basename, sentence)).strip('\n'))


    def _frag_restart(self):
        """
        To restart FragPELE fragment libraries simulation."
        """
        self._extract_working_directory()
        with open(self.parameters.input, "r") as conf_file:
            final_bonds = []
            for line, path in zip(conf_file.readlines(), self.parameters.working_dir):
                if len(glob.glob(os.path.join(path, "*_top.pdb"))) == 0:
                    if os.path.exists(path):
                        shutil.rmtree(path)
                    final_bonds.append(line)
            with open(self.parameters.input, "w") as conf_file:
                for line in final_bonds:
                    conf_file.write(line + "\n")