import glob
import os

from pele_platform.PPI.cluster import cluster_best_structures
from pele_platform.PPI.simulation_launcher import launch_simulation
from pele_platform.PPI.preparation import prepare_structure
from pele_platform.Utilities.Helpers.helpers import cd, get_next_peledir, get_latest_peledir
import pele_platform.Utilities.Parameters.parameters as pv
import pele_platform.Adaptive.simulation as si
from pele_platform.Adaptive import parametrizer


def run_ppi(parsed_yaml: dict) -> (pv.ParametersBuilder, pv.ParametersBuilder):
    # Let user choose working folder
    original_dir = os.path.abspath(os.getcwd())
    working_folder = os.path.abspath("{}_Pele".format(parsed_yaml.residue))
    if not parsed_yaml.folder:
        working_folder = get_next_peledir(working_folder) if not parsed_yaml.adaptive_restart else get_latest_peledir(
            working_folder)
    else:
        working_folder = os.path.abspath(parsed_yaml.folder)

    # Set main folder
    parsed_yaml.folder = os.path.join(working_folder, "1_interface_exploration")

    # get arguments from input.yaml
    n_waters = parsed_yaml.n_waters
    parsed_yaml.n_waters = None
    protein_file = parsed_yaml.system
    chain = parsed_yaml.protein
    ligand_pdb = parsed_yaml.ligand_pdb

    # Parametrize hetero molecules before merging PDBs, if using peleffy. Otherwise they will go through Plop in
    # Adaptive.simulation.
    if parsed_yaml.use_peleffy:
        templates, rotamers, to_skip = parametrize_hetero_ppi(parsed_yaml)

        parsed_yaml.template = parsed_yaml.template.extend(templates) if parsed_yaml.template else templates
        parsed_yaml.rotamers = parsed_yaml.rotamers.extend(rotamers) if parsed_yaml.rotamers else rotamers
        parsed_yaml.skip_ligand_prep = parsed_yaml.skip_ligand_prep.extend(to_skip) if parsed_yaml.skip_ligand_prep else to_skip

    # no waters in the first simulation
    parsed_yaml.water_arg = None
    parsed_yaml.use_peleffy = parsed_yaml.use_peleffy if parsed_yaml.use_peleffy is not None else False

    # remove chains except for "protein" flag
    protein_file = prepare_structure(protein_file, ligand_pdb, chain, True, peleffy=parsed_yaml.use_peleffy)
    parsed_yaml.system = protein_file

    # start simulation 1 - induced fit
    parsed_yaml.induced_fit_exhaustive = True
    simulation1 = si.run_adaptive(parsed_yaml)
    simulation1_path = os.path.join(simulation1.pele_dir, simulation1.output)

    # cluster best structures
    if not parsed_yaml.debug:
        with cd(simulation1_path):
            cluster_best_structures("5", n_components=simulation1.n_components,
                                    residue=simulation1.residue, topology=simulation1.topology,
                                    directory=working_folder, logger=simulation1.logger)

    # adjust original input.yaml
    if not parsed_yaml.skip_refinement:
        parsed_yaml.system = os.path.join(working_folder, "refinement_input/*.pdb")
        parsed_yaml.folder = os.path.join(working_folder, "2_refinement_simulation")
        parsed_yaml.induced_fit_exhaustive = None
        parsed_yaml.ppi = None
        parsed_yaml.poses = None
        parsed_yaml.rescoring = True
        del parsed_yaml.water_arg

        # Set waters only if specified by user
        if n_waters != 0:
            parsed_yaml.waters = "all_waters"
            parsed_yaml.n_waters = n_waters
        else:
            parsed_yaml.waters = None
            parsed_yaml.n_waters = n_waters
        parsed_yaml.adaptive_restart = False
        if not parsed_yaml.test:
            parsed_yaml.iterations = 1
            parsed_yaml.steps = 100
        parsed_yaml.box_center = simulation1.box_center
        parsed_yaml.box_radius = 100  # We should have a look at how to set no box but at the moment super big

        # start simulation 2 - minimisation
        with cd(original_dir):
            if not parsed_yaml.debug:
                simulation2 = launch_simulation(parsed_yaml)
            else:
                simulation2 = None
    else:
        simulation2 = None
    return simulation1, simulation2


def parametrize_hetero_ppi(parsed_yaml):
    """
    Parametrizes hetero molecules in ligand and protein files before the two get merged, then retrieve created files
    and ligands to skip.

    Parameters
    ------------
    parsed_yaml : dict
        Parsed input YAML in a form of a dictionary.

    Returns
    ---------
    template_files : List[str]
        List of paths to created template files.
    rotamer_files : List[str]
        List of paths to created rotamer files.
    to_skip : List[str]
        List of residues to skip - this ensures Parametrized will not try to process them later in Adaptive.simulation.
    """
    # Fix mismatch between YamlParser and Parameters attributes
    parsed_yaml.external_template = parsed_yaml.template
    parsed_yaml.external_rotamers = parsed_yaml.ext_rotamers
    parsed_yaml.pele_dir = "."
    parsed_yaml.as_datalocal = False

    # Parametrize ligand and protein
    obj = parametrizer.Parametrizer.from_parameters(parsed_yaml)
    obj.parametrize_ligands_from(parsed_yaml.ligand_pdb)
    obj.parametrize_ligands_from(parsed_yaml.system)

    # Retrieve created files
    template_files = glob.glob(os.path.join(os.getcwd(), "*z"))
    rotamer_files = glob.glob(os.path.join(os.getcwd(), "*.rot.assign"))

    # Figure out which residues should be skipped
    to_skip = [os.path.basename(file).rstrip("z") for file in template_files]
    to_skip.extend([os.path.basename(file).rstrip("rot.assign").lower() for file in rotamer_files])

    # Remove duplicates
    to_skip = list(set(to_skip))

    return template_files, rotamer_files, to_skip
