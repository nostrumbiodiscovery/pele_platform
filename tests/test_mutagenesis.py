import glob
import os
import shutil

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

    # List of mutated PDBs expected in each job
    expected_pdbs = [('T454A_processed.pdb', 'T454D_processed.pdb'),
                     ('T454E_processed.pdb', )]

    for i, job in enumerate(all_jobs):
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
        job_dir = os.path.join(job.pele_dir, "input", "*.pdb")
        job_input = [os.path.basename(file) for file in glob.glob(job_dir)]
        for expected_pdb in expected_pdbs[i]:
            assert expected_pdb in job_input


def test_mutagenesis_restart():
    yaml = os.path.join(test_path, "restart_saturated_mutagenesis.yaml")
    pele_dir = "restart_ANL_Pele"
    restart_folder = os.path.join(test_path, "directory_to_restart")

    if os.path.exists(pele_dir):
        shutil.rmtree(pele_dir)
    shutil.copytree(restart_folder, pele_dir)

    all_jobs = main.run_platform(yaml)
    last_job = all_jobs[-1]
    mutation_folder = os.path.splitext(os.path.basename(last_job.input[0]))[0]
    restart_trajectories = glob.glob(os.path.join(last_job.pele_dir, last_job.output, mutation_folder, "traj*"))
    assert len(restart_trajectories) == 2
