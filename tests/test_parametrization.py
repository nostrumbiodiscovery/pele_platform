import os
import pytest

from pele_platform.Adaptive import parametrization
from pele_platform.constants import constants

test_path = os.path.join(constants.DIR, "Examples", "constraints")


@pytest.mark.parametrize(
    ("pdb", "ligands"),
    [
        (os.path.join(test_path, "1zop.pdb"), ["MG", "LIG"]),
        (os.path.join(test_path, "4zu9_prep.pdb"), ["MET", "MG", "SO4", "LIG"]),
    ],
)
def test_ligand_extraction(pdb, ligands):
    """
    Checks hetatm ligand extraction from PDB file in Parametrization class.

    Parameters
    -----------
    pdb : str
        Path to PDB file.
    ligands : List[str]
        List of expected residue names to be extracted from the PDB.
    """
    extracted_ligands = parametrization.Parametrization.extract_ligands(pdb, gridres=10)
    assert len(extracted_ligands) == len(ligands)


@pytest.mark.xfail(reason="Waiting for Laura to fix spaces in file names.")
def test_generate_ligand_parameters(parametrize):
    """
    Runs the whole parametrization flow and checks if all expected files were created.

    Parameters
    -----------
    parametrize : Parametrization
        Object created in fixture.
    """
    expected_templates = ["mgz", "ligz"]
    expected_rotamers = ["MG.rot.assign", "LIG.rot.assign"]
    parametrize.generate_ligand_parameters()

    for template in expected_templates:
        assert os.path.exists(template)

    for rotamer in expected_rotamers:
        assert os.path.exists(rotamer)


@pytest.mark.parametrize(
    ("solvent", "forcefield", "error"),
    [
        ("OBC2", "OPLS2005", False),
        ("vdgbnp", "OPLS2005", False),
        ("OBC2", "openff-1.2.0", False),
        ("vdgbnp", "openff-1.3.0", True),
    ],
)
def test_check_solvent(solvent, forcefield, error):
    """
    Checks if certain combinations of forcefield and solvent raise error when necessary.

    Parameters
    -----------
    solvent : str
        User-defined solvent, one of parametrization.forcefields.
    forcefield : str
        User-defined solvent, one of parametrization.solvents
    error : bool
        Set to True if particular combination is supposed to raise en error.
    """
    if error:
        with pytest.raises(ValueError):
            parametrization.Parametrization._check_solvent(
                solvent=solvent, forcefield=forcefield
            )
    else:
        parametrization.Parametrization._check_solvent(
            solvent=solvent, forcefield=forcefield
        )


@pytest.fixture
def parametrize():
    """
    Creates Parametrization class to use in other tests.
    """
    obj = parametrization.Parametrization(
        pdb_file=os.path.join(test_path, "1zop.pdb"),
        forcefield="openff-1.3.0",
        solvent="OBC", as_datalocal=False)
    return obj
