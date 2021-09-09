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
import frag_pele.main as frag
from pele_platform.analysis.analysis import Analysis
from pele_platform.Utilities.Helpers import water
from pele_platform.context import context


class FragRunner(object):
    def __init__(self):
        """
        Initializes FragRunner class.
        """
        context.parameters_builder.build_frag_variables()

    def run(self):
        self._set_test_variables()
        self._prepare_control_file()
        self._launch()
        if not context.parameters.debug:
            self._analysis()
        return context.parameters

    def _launch(self):
        params = context.parameters
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
        params = context.parameters

        tmp_dir = tempfile.mkdtemp()
        tmp_control_file = os.path.join(tmp_dir, os.path.basename(params.control_file))
        shutil.copy(params.control_file, tmp_control_file)
        adaptive = ad.SimulationBuilder("", tmp_control_file, params)

        # Fill to time because we have flags inside flags
        formatted_params = adaptive.format_parameters()
        adaptive.fill_pele_template(formatted_params, params.water_object)
        adaptive.fill_pele_template(formatted_params, params.water_object)
        context.parameters.control_file = tmp_control_file

        return context.parameters.control_file

    def _prepare_parameters(self):
        context.parameters.spython = cs.SCHRODINGER  # Commented to use Frag 2.2.1 instead of Frag 3.0.0

    def _set_test_variables(self):
        if context.parameters.test:
            context.parameters.gr_steps = 1
            context.parameters.frag_steps = 2
            context.parameters.frag_eq_steps = 1
            context.parameters.temperature = 100000
            context.parameters.anm_freq = 0
            context.parameters.minimizatoon = 0
            context.parameters.sidechain_freq = 0
            context.parameters.water_freq = 0
            context.parameters.cpus = 4
        else:
            pass

    def _run(self):
        params = context.parameters
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

        context.parameters.input = "input.conf"

        params = context.parameters

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

        params = context.parameters
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
        if context.parameters.only_analysis:
            self._extract_working_directory()
            for path in context.parameters.working_dir:
                if not os.path.exists(path):
                    if not os.path.exists(context.parameters.folder):
                        raise ce.MissingWorkingDir("Missing working_folder parameter. Please set the "
                                                   "path of the trajectories using the flag working_"
                                                   "folder in your input.yaml")
            self._prepare_parameters()
            context.parameters.water_ids_to_track = water.water_ids_from_conf(context.parameters.control_file)

        if context.yaml_parser.analysis_to_point:
            context.parameters.analysis_to_point = context.yaml_parser.analysis_to_point
            ana.main(
                path=context.parameters.folder,
                atom_coords=context.parameters.analysis_to_point,
                pattern=os.path.basename(context.parameters.system),
            )
        for path in context.parameters.working_dir:
            simulation_output = os.path.join(path, 'sampling_result')
            analysis_folder = os.path.join(path, "results")

            analysis = Analysis(
                resname="GRW",
                chain=context.parameters.chain,
                simulation_output=simulation_output,
                be_column=context.parameters.be_column,
                limit_column=context.parameters.limit_column,
                traj=context.parameters.traj_name,
                report=context.parameters.report_name,
                skip_initial_structures=not context.parameters.test,
                kde=context.parameters.kde,
                kde_structs=context.parameters.kde_structs,
                topology=context.parameters.topology,
                cpus=context.parameters.cpus,
                water_ids_to_track=context.parameters.water_ids_to_track,
            )
            analysis.generate(
                analysis_folder,
                clustering_type=context.parameters.clustering_method.lower(),
                bandwidth=context.parameters.bandwidth,
                analysis_nclust=context.parameters.analysis_nclust,
                max_top_clusters=context.parameters.max_top_clusters,
                top_clusters_criterion=context.parameters.top_clusters_criterion,
                min_population=context.parameters.min_population,
                representatives_criterion=context.parameters.cluster_representatives_criterion,
            )

    def _clean_up(self, fragment_files):
        for file in fragment_files:
            if os.path.isfile(file):
                os.remove(file)

    def _extract_working_directory(self):
        params = context.parameters
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
                if context.parameters.folder:
                    params.working_dir.append(
                        os.path.join(context.parameters.folder, "{}_{}".format(pdb_basename, sentence)).strip('\n'))
                else:
                    params.working_dir.append(
                        os.path.join(current_path, "{}_{}".format(pdb_basename, sentence)).strip('\n'))

    def _frag_restart(self):
        self._extract_working_directory()
        with open(context.parameters.input, "r") as conf_file:
            final_bonds = []
            for line, path in zip(conf_file.readlines(), context.parameters.working_dir):
                if len(glob.glob(os.path.join(path, "*_top.pdb"))) == 0:
                    if os.path.exists(path):
                        shutil.rmtree(path)
                    final_bonds.append(line)
            with open(context.parameters.input, "w") as conf_file:
                for line in final_bonds:
                    conf_file.write(line + "\n")
