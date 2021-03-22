import os
import logging
import numpy as np
import sys
import warnings
import PPP.global_variables as gv
from Bio.PDB import PDBParser
import pele_platform.Errors.custom_errors as cs
from multiprocessing import Pool
from functools import partial


__all__ = ["get_suffix", "backup_logger"]


def silentremove(*args, **kwargs):
    for files in args:
        for filename in files:
            try:
                os.remove(filename)
            except OSError:
                pass


def create_dir(base_dir, extension=None):
    """
        Class Method to manage
        directory creation only if that
        ones doesn't exist

        Location:
            base_dir+extension
            or base_dir if extension is None
    """
    if extension:
        path = os.path.join(base_dir, extension)
        if os.path.isdir(path):
            warnings.warn("Directory {} already exists.".format(path), RuntimeWarning)
        else:
            os.makedirs(path)
    else:
        if os.path.isdir(base_dir):
            warnings.warn("Directory {} already exists.".format(base_dir), RuntimeWarning)
        else:
            os.makedirs(base_dir)


class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)


def is_repeated(pele_dir):
    """
    Given a PELE directory it will return a new directory with a new
    suffix. The suffix is chosen with the following criterion:

     - In case that NOL_Pele folder already exists in the working directory,
       it will return NOL_Pele_1.
     - In case that NOL_Pele_1 folder already exists in the working directory,
       it will return NOL_Pele_2.

    .. todo ::
        * The name of this function is misleading.

    Parameters
    ----------
    pele_dir : str
        The candidate name for the PELE directory which will only be
        modified following the criterion above if it already exists

    Returns
    -------
    pele_dir : str
        The new PELE directory that does not match with any other directory
        previously created
    """

    original_dir = None
    split_dir = pele_dir.split("_")
    for chunk in split_dir:
        if chunk != "Pele":
            if original_dir:
                original_dir = "{}_{}".format(original_dir, chunk)
            else:
                original_dir = chunk
        else:
            break
    if split_dir[-1].isdigit():
        i = split_dir[-1]
        i = int(i) + 1
    else:
        i = 1
    if os.path.isdir(pele_dir):
                new_pele_dir = "{}_Pele_{}".format(original_dir, i)
                new_pele_dir = is_repeated(new_pele_dir)
                return new_pele_dir
    else:
                return pele_dir


def is_last(pele_dir):
    """
    Given a PELE directory it will return the name of the directory that
    looks newer. It employs the following criterion:

     - In case that NOL_Pele, NOL_Pele_1 and NOL_Pele_2 folders already
       exist in the working directory, it will choose NOL_Pele_2 since it
       has the highest suffix index.
     - In case no directory named NOL_Pele is found, the original name
       will be employed.

    .. todo ::
        * The name of this function is misleading.

    Parameters
    ----------
    pele_dir : str
        The original name for the PELE directory whose newest directory
        wants to be obtained

    Returns
    -------
    pele_dir : str
        The newest PELE directory that has been found in the working
        directory according to the original name that is supplied
    """

    original_dir = None
    split_dir = pele_dir.split("_")
    for chunk in split_dir:
                if chunk != "Pele":
                        if original_dir:
                                original_dir = "{}_{}".format(original_dir, chunk)
                        else:
                                original_dir = chunk
                else:
                        break
    if split_dir[-1].isdigit():
        i = split_dir[-1]
        i = int(i) + 1
    else:
                i = 1

    if os.path.isdir(pele_dir):
            new_pele_dir = "{}_Pele_{}".format(original_dir, i)
            if not os.path.isdir(new_pele_dir):
                return pele_dir
            else:
                            new_pele_dir = is_last(new_pele_dir)
                            return new_pele_dir
    else:
        return pele_dir


def retrieve_atom_info(atom, pdb):
    """
    Parse pdb and return atom name
    chain and residue number
    """
    with open(pdb, "r") as f:
        for line in f:
            try:
                if not isinstance(atom, int) and not atom.isdigit():
                    try:
                        chain, resnum, atomname = atom.split(":")
                    except ValueError:
                        raise ValueError(f"Check atom distance entrance {atom}. Should be like this: 'A:220:OD1'")
                    if line[21].strip() == chain and line[22:26].strip() == resnum and line[12:16].strip() == atomname:
                        atomname = line[12:16]
                        return chain + ":" + resnum + ":" + atomname.replace(" ", "_")
                else:
                    if line[6:11].strip() == str(atom):
                        chain = line[21].strip() 
                        resnum = line[22:26].strip()
                        atomname = line[12:16]
                        return chain + ":" + resnum + ":" + atomname.replace(" ", "_")
            except IndexError:
                pass
        sys.exit(f"Check the atoms {atom} given to calculate the distance metric.")


def retrieve_all_waters(pdb, exclude=False):
    with open(pdb, 'r') as f:
        waters = list(set(["{}:{}".format(line[21:22], line[22:26].strip()) for line in f if line and "HOH" in line]))
    if exclude:
        waters = [water for water in waters if water not in exclude]
    return waters


def retrieve_constraints_for_pele(constraints, pdb):
    CONSTR_ATOM_POINT = '{{ "type": "constrainAtomToPosition", "springConstant": {}, "equilibriumDistance": 0.0, "constrainThisAtom": "{}:{}:{}" }},'
    CONSTR_ATOM_ATOM = '{{"type": "constrainAtomsDistance", "springConstant": {}, "equilibriumDistance": {}, "constrainThisAtom":  "{}:{}:{}", "toThisOtherAtom": "{}:{}:{}"}},'
    final_constraints = []
    for constraint in constraints:
        # Atom to point constraint: 2.2-A:123:2 or 2.2-1986
        if len(constraint.split("-")) == 2:
            spring_constant, atom_info = constraint.split("-")
            chain, residue, atom_name = retrieve_atom_info(atom_info, pdb).split(":")
            constraint = CONSTR_ATOM_POINT.format(spring_constant, chain, residue, atom_name)
        # Atom to atom constraint: 2.2-2.75-A:123:2-B:2:7 or 2.2-2.74-1985-1962
        elif len(constraint.split("-")) == 4:
            spring_constant, eq_distance, atom1_info, atom2_info = constraint.split("-")
            chain1, residue1, atom_name1 = retrieve_atom_info(atom1_info, pdb).split(":")
            chain2, residue2, atom_name2 = retrieve_atom_info(atom2_info, pdb).split(":")
            constraint =  CONSTR_ATOM_ATOM.format(spring_constant, eq_distance, chain1, residue1, atom_name1, chain2, residue2, atom_name2)
        final_constraints.append(constraint)
    return final_constraints


def retrieve_box(structure, residue_1, residue_2, weights=[0.5, 0.5]):
    # get center of interface (if PPI)
    coords1 = get_coords_from_residue(structure, residue_1)
    coords2 = get_coords_from_residue(structure, residue_2)
    
    box_center = np.average([coords1, coords2], axis=0, weights=weights)
    box_radius = abs(np.linalg.norm(coords1-coords2))/2 + 4 #Sum 4 to give more space
    return list(box_center), box_radius


def get_coords_from_residue(structure, original_residue):
    parser = PDBParser()
    structure = parser.get_structure('protein', structure)
    chain, res_number, atom_name = original_residue.split(":")
    try:
        res_number = int(res_number)
    except ValueError:
        raise cs.WrongAtomStringFormat(f"The specified atom is wrong '{original_residue}'. \
Should be 'chain:resnumber:atomname'")
    for residue in structure.get_residues():
        if residue.id[1] == res_number:
            for atom in residue.get_atoms():
                if atom.name == atom_name:
                    COI = np.array(list(atom.get_vector()))
                    return COI
    raise cs.WrongAtomSpecified(f"Atom {original_residue} could not be found in structure")


def backup_logger(logger, message):
    if not logger:
        logger = logging.getLogger('logger')
        logger.setLevel(logging.INFO)
        logger.info(message)
    else:
        logger.info(message)


def find_nonstd_residue(pdb):
    with open(pdb, "r") as f:
        resnames = list(set([line[17:20] for line in f \
    if line.startswith("ATOM") and line[17:20] not in gv.supported_aminoacids]))
        return resnames


def parallelize(func, iterable, n_workers, **kwargs):
    pool = Pool(n_workers)
    f = partial(func, **kwargs)
    return pool.map(f, iterable)
    pool.close()
    pool.join()


def is_rdkit():
    try:
        import rdkit
        from rdkit import Chem
        return True
    except:
        raise ModuleNotFoundError("Please install rdkit with the following command: conda install -c conda-forge rdkit")


def get_suffix(filename, separator='_'):
    """
    Given a filename, it returns its corresponding suffix.

    Parameters
    ----------
    filename : str
        The filename path
    separator : str
        The pattern that is used to separate the name root from the suffix.
        Default is '_'

    Returns
    -------
    suffix : str
        The suffix for the supplied filename
    """
    name = os.path.basename(filename)
    suffix = name.split('_')[-1]

    return suffix