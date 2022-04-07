from copy import deepcopy
from dataclasses import dataclass, field
import glob
import os
from typing import List

from pele_platform.Utilities.Parameters import parameters
from pele_platform.Utilities.Helpers import helpers, yaml_parser
from pele_platform.Adaptive import simulation
from pele_platform.enzyme_engineering.saturated_mutagenesis import SaturatedMutagenesis

from satumut.simulation import SimulationRunner
from satumut.helper import neighbourresidues
from satumut.mutate_pdb import generate_mutations
from satumut.PELE_Plurizymer_analysis import PELE_PLURIZYMER_add_metrics
from satumut.pele_files import CreateYamlFiles

@dataclass
class Plurizymer(SaturatedMutagenesis):
    """
    Interface to run plurizymer simulation from another repository.

    Parameters
    ----------
    env : yaml_parser.YamlParser
        Arguments provided by the user in input.yaml.
    already_computed : List[str]
        Initially empty list of already computed systems.
    all jobs: List[EnviroBuilder]
        Initially empty list of all completed jobs.
    original_dir : str
        Directory from which the job is launched.
    start : int
        Index to enumerate subset folders, if restarting adaptive,
        otherwise default = 1
    subset_folder : str
    See Also folder name, default = "Subset_"
    """

    env: yaml_parser.YamlParser
    already_computed: List = field(default_factory=list)
    all_jobs: List = field(default_factory=list)
    original_dir: str = os.path.abspath(os.getcwd())
    start: int = 1
    subset_folder: str = "Subset_"

    def __post_init__(self):
        """
            Parse the input parameters
        """
        builder = parameters.ParametersBuilder()
        self.params = builder.build_plurizymer_variables(self.env)

    def run(self):
        """
        Runs the simulation for all inputs.

        Returns
        -------
        all_jobs : list
            A list of job parameters (EnviroBuilder objects) for each
            simulation subset
        """
        self.set_package_params()
        self.check_cpus()
        self.set_working_folder()
        self.params.folder = f"{self.working_folder}_mut"
        simulation_satumut = SimulationRunner(self.params.system, self.params.folder)
        input_ = simulation_satumut.side_function()
        # use this object to keep track of the folder where the simulations
        # will be stored
        satumut_helper = CreateYamlFiles([], "", "", cpus=self.params.cpus,
                                         single=self.params.plurizymer_single_mutation,
                                         turn=self.params.plurizymer_turn)
        position = neighbourresidues(input_, self.params.plurizymer_atom,
                                     self.params.satumut_radius_neighbors,
                                     self.params.satumut_fixed_residues)
        if not self.params.adaptive_restart:
            generate_mutations(input_, position, hydrogen=self.params.satumut_hydrogens,
                               pdb_dir=self.params.satumut_pdb_dir,
                               single=self.params.plurizymer_single_mutation,
                               turn=self.params.plurizymer_turn)
        self.all_mutations = [os.path.abspath(x) for x in glob.glob(os.path.join(self.params.satumut_pdb_dir, "*.pdb"))]
        self.check_metric_distance_atoms(input_)

        self.set_simulation_folder(satumut_helper)
        self.restart_checker()
        self.split_into_subsets()


        for idx, subset in enumerate(self.mutation_subsets, self.start):
            self.env.input = subset
            self.env.cpus = self.env.cpus_per_mutation * len(subset) + 1
            self.env.folder = os.path.join(
                self.working_folder, "{}{}".format(self.subset_folder, idx)
            )

            with helpers.cd(self.original_dir):
                job = simulation.run_adaptive(self.env)
                self.postprocessing(job)
                self.all_jobs.append(deepcopy(job))
                self.logger(job)

        PELE_PLURIZYMER_add_metrics(self.params.folder,
                self.params.plurizymer_single_mutation,
                self.params.plurizymer_atom,
                self.params.xtc, self.params.cpus)

        os.chdir(self.original_dir)
        return self.all_jobs
