import os
import shutil

from pele_platform.Utilities.Helpers.ligand_conformations import LigandConformations
from pele_platform.constants import constants as cs
from pele_platform import main
from .utils import check_file

test_path = os.path.join(cs.DIR, "Examples", "ligand_conformations")


expected_lines = [
    '"ConformationPerturbation":{',
    '"overlapFactor": 0.66',
    '"conformationPerturbationFrequency": 2',
]


def test_class():
    """
    Tests class that generates conformation library using 3 ligand PDB files as input.
    """
    expected_output = "DataLocal/Conformations/IK1.conformation"
    system_file = os.path.join(cs.DIR, "Examples", "out_in", "1_3ZON_complex.pdb")

    generator = LigandConformations(
        path=test_path,
        system=system_file,
        resname="IK1",
        forcefield="OPLS2005",
    )
    generator.generate()

    assert os.path.exists(expected_output)
    shutil.rmtree(expected_output.split("/")[0])


def test_production():
    """
    Tests the full flow in debug mode and checks for expected values in pele.conf as well as generated libraries.
    """
    yaml = os.path.join(test_path, "input.yaml")
    job = main.run_platform_from_yaml(yaml)

    errors = check_file(job.pele_dir, "pele.conf", expected_lines, [])
    assert not errors

    assert os.path.exists(os.path.join(job.pele_dir, "DataLocal", "Conformations", f"{job.residue}.conformation"))
