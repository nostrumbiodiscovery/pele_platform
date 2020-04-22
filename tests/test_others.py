import os
import glob
import pele_platform.constants.constants as cs
import pele_platform.constants.pele_params as pp
import pele_platform.main as main
from . import test_adaptive as tk
import pytest
import pele_platform.Utilities.Helpers.protein_wizard as pp
import pele_platform.Frag.checker as ch
import pele_platform.Errors.custom_errors as ce

test_path = os.path.join(cs.DIR, "Examples")
EXTERNAL_CONSTR_ARGS = os.path.join(test_path, "constraints/input_external_constraints.yaml")
PPP_CONSTR_ARGS = os.path.join(test_path, "constraints/input_ppp.yaml")


EXT_CONSTR = [
    '{ "type": "constrainAtomToPosition", "springConstant": 10, "equilibriumDistance": 0.0, "constrainThisAtom": "A:1:_H__" },',
    '{ "type": "constrainAtomToPosition", "springConstant": 5, "equilibriumDistance": 0.0, "constrainThisAtom": "A:1:_H__" },',
    '{"type": "constrainAtomsDistance", "springConstant": 50, "equilibriumDistance": 2.34, "constrainThisAtom":  "A:1:_H__", "toThisOtherAtom": "L:1:_C21"},',
    '{"type": "constrainAtomsDistance", "springConstant": 50, "equilibriumDistance": 2.34, "constrainThisAtom":  "A:1:_H__", "toThisOtherAtom": "L:1:_C21"}'
]

PPP_CONSTR = [
    '"constrainThisAtom": "B:207:_CA_" }',
    '"constrainThisAtom": "A:1:_CA_" }',
    '"constrainThisAtom": "B:247:_CA_" }'
]

def test_external_constraints(ext_args=EXTERNAL_CONSTR_ARGS):
    errors = []
    job = main.run_platform(ext_args)
    errors = tk.check_file(job.pele_dir, "pele.conf", EXT_CONSTR, errors)
    assert not errors


def test_ppp_constraints(ext_args=PPP_CONSTR_ARGS):
    errors = []
    job = main.run_platform(ext_args)
    errors = tk.check_file(job.pele_dir, "pele.conf", PPP_CONSTR, errors)
    assert not errors


def test_checker():
    yaml = os.path.join(test_path, "checker/input.yaml")
    try:
        job = main.run_platform(yaml)
    except KeyError:
        assert KeyError
        return
    assert False


def test_proteinwizard():
    complex_correct = os.path.join(test_path, "preparation/6qmk_correct.pdb")
    complex_repeated = os.path.join(test_path, "preparation/6qmk_repeated.pdb")
    yaml = os.path.join(test_path, "preparation/input.yaml")

    correct_output = pp.prep_complex(complex_correct, yaml, prep_output=complex_correct, debug=True)
    repeated_output = pp.prep_complex(complex_repeated, yaml, prep_output=complex_repeated, debug=True)
    
    with open(correct_output, "r") as f:
        lines = f.readlines()
        for line in lines:
            correct_lig_chains = [line[21:22].strip() for line in lines if line[17:20] == "J8H"]
            correct_lig_atomnames = [line[12:16].strip() for line in lines if line [17:20] == "J8H"]

    with open(repeated_output, "r") as f:
        lines = f.readlines()
        for line in lines:
            repeated_lig_chains = [line[21:22].strip() for line in lines if line[17:20] == "J8H"]
            repeated_lig_atomnames = [line[12:16].strip() for line in lines if line [17:20] == "J8H"]

    assert all(elem == "Z" for elem in repeated_lig_chains) and all(elem == "Z" for elem in correct_lig_chains) # make sure all lig chains are Z
    assert len(set(correct_lig_atomnames)) == len(correct_lig_atomnames) # check for repetition
    assert len(set(repeated_lig_atomnames)) == len(repeated_lig_atomnames)


SDF_LIMIT_ATOMS = os.path.join(test_path, "frag/checker/receptor.sdf")
def test_checker_nlimits(sdf=SDF_LIMIT_ATOMS):
    try:
        ch.check_limit_number_atoms(sdf, 100)
    except ce.LigandSizeExceed:
        assert True
    else:
        assert False
    

PDB_CHECKER_SUBSEARCH = os.path.join(test_path, "frag/checker/4RFM_series.sdf")
CORE_CHECKER_SUBSEARCH = os.path.join(test_path, "frag/checker/4RFM_proc.pdb")
def test_checker_subsearch(ligand=PDB_CHECKER_SUBSEARCH, core=CORE_CHECKER_SUBSEARCH):
    from rdkit import Chem
    mol = next(Chem.SDMolSupplier(ligand))
    core = Chem.rdmolfiles.MolFromPDBFile(core)
    atoms_in_common = mol.GetSubstructMatches(core)[0]
    atoms_in_common_after = ch.chec_substructure_match(core, mol, atoms_in_common)
    # Exchange nitrogen due to wrong previous result
    assert atoms_in_common != atoms_in_common_after
    assert atoms_in_common_after[atoms_in_common.index(13)] == 12
    assert atoms_in_common_after[atoms_in_common.index(12)] == 13
    
