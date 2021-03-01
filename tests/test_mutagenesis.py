import glob
import os

from pele_platform.constants import constants
from pele_platform import main
from tests import test_adaptive_defaults

test_path = os.path.join(constants.DIR, "Examples/enzyme_engineering")


def test_mutagenesis_production():
    """
    Tests end-to-end saturated mutagenesis run on 5 CPUs to make sure we get the right output and pele.conf looks
    as expected (i.e. has induced_fit_exhaustive defaults).
    TODO: This is a very primitive implementation, needs a lot more testing!
    """
    yaml = os.path.join(test_path, "saturated_mutagenesis.yaml")
    all_jobs = main.run_platform(yaml)

    for job in all_jobs:
        # Output files exist
        output_files = glob.glob(os.path.join(job.pele_dir, job.output, "*/traj*.pdb"))
        assert output_files

        # Check if the default induced fit exhaustive lines are present in pele configuration file
        induced_fit_lines = test_adaptive_defaults.INDUCE_FIT_PELE
        errors = test_adaptive_defaults.check_file(job.pele_dir, "pele.conf", induced_fit_lines, [])
        assert not errors

        # Make sure all subset directories have correct names
        assert "Subset" in os.path.basename(job.pele_dir)

        # Analysis should not run
        results_folder = os.path.exists(os.path.join(job.pele_dir, "results"))
        assert not results_folder

    # Check if the jobs were properly split between subsets based on available CPUs
    last_job_dir = os.path.join(all_jobs[-1].pele_dir, "input", "*.pdb")
    last_job_input = [os.path.basename(file) for file in glob.glob(last_job_dir)]
    assert ('T454W_processed.pdb' and 'original_processed.pdb') in last_job_input
