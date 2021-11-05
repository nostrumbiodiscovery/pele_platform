import os
import glob

from pele_platform.constants import constants
from pele_platform import main
from tests import test_adaptive_defaults

test_path = os.path.join(constants.DIR, "Examples/enzyme_engineering")


def test_plurizymer_production():
    """
    Tests end-to-end plurizymer run on make sure we get the right output and pele.conf looks
    as expected (i.e. has induced_fit_exhaustive defaults).
    """
    yaml = os.path.join(test_path, "plurizymer.yaml")
    all_jobs = main.run_platform_from_yaml(yaml)

    # List of mutated PDBs expected in each job
    expected_pdbs = [
        ("N209S_processed.pdb", "original_processed.pdb", "T212S_processed.pdb"),
        ("L215S_processed.pdb",),
    ]

    for i, job in enumerate(all_jobs):
        # Output files exist
        output_files = glob.glob(os.path.join(job.pele_dir, job.output, "*/traj*.pdb"))
        assert output_files

        # Check if all files from iterations folder were postprocessed and moved to their respective mutation folders
        assert not glob.glob(os.path.join(job.pele_dir, job.output, "0/traj*pdb"))

        # Check if the default induced fit exhaustive lines are present in pele configuration file
        induced_fit_lines = test_adaptive_defaults.INDUCE_FIT_PELE
        errors = test_adaptive_defaults.check_file(
            job.pele_dir, "pele.conf", induced_fit_lines, []
        )
        assert not errors

        # Make sure all subset directories have correct names
        assert "Subset" in os.path.basename(job.pele_dir)

        # Analysis should not run
        results_folder = os.path.exists(os.path.join(job.pele_dir, "results"))
        assert not results_folder

        # Check if the jobs were properly split between subsets based on available CPUs
        job_dir = os.path.join(job.pele_dir, "input", "*.pdb")
        job_input = [os.path.basename(file) for file in glob.glob(job_dir)]
        for expected_pdb in expected_pdbs[i]:
            assert expected_pdb in job_input
