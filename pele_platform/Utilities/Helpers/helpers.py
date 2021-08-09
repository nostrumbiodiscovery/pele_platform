import os
import logging
import numpy as np
import re
import subprocess
import shutil
import sys
import warnings
import PPP.global_variables as gv
from Bio.PDB import PDBParser
import pele_platform.Errors.custom_errors as cs
from multiprocessing import Pool
from functools import partial
from packaging import version
from pele_platform.Errors.custom_errors import ResidueNotFound, PELENotFound


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
    It creates a directory only if that one doesn't exist.

    Parameters
    ----------
    base_dir : str
        The base directory to remove
    extension : str
        The specific extension to remove, if any. Default is None
    """
    if extension:
        path = os.path.join(base_dir, extension)
        if os.path.isdir(path):
            warnings.warn("Directory {} already exists.".format(path), RuntimeWarning)
        else:
            os.makedirs(path)
    else:
        if os.path.isdir(base_dir):
            warnings.warn(
                "Directory {} already exists.".format(base_dir), RuntimeWarning
            )
        else:
            os.makedirs(base_dir)


class cd:
    """Context manager for changing the current working directory."""

    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)


def get_directory_new_index(pele_dir):
    """
    It gets the new index of the next PELE directory to be created,
    for example given LIG_Pele_11 will return new_index = 12 and
    original_dir = "LIG".

    Parameters
    ----------
    pele_dir : str
        The PELE directory, e.g. LIG_Pele_1

    Returns
    -------
    new_index : int
        The index corresponding to the next PELE directory
    old_index : int
        The index corresponding to the current PELE directory
    original_dir : str
        The original PELE directory (usually, it matches with the
        residue name)
    """
    folder_name = os.path.basename(pele_dir)
    split_name = folder_name.split("_")

    if len(split_name) < 2 or len(split_name) > 3 or split_name[1] != "Pele":
        raise ValueError(
            "Invalid PELE directory {}. ".format(folder_name) + "Its format is unknown"
        )

    original_dir = split_name[0]

    if split_name[-1].isdigit():
        new_index = split_name[-1]
        new_index = int(new_index) + 1
    else:
        new_index = 1

    return new_index, new_index - 1, original_dir


def get_next_peledir(pele_dir):
    """
    Given a PELE directory it will return a new directory with a new
    suffix. The suffix is chosen with the following criterion:

     - In case that NOL_Pele folder already exists in the working directory,
       it will return NOL_Pele_1.
     - In case that NOL_Pele_1 folder already exists in the working directory,
       it will return NOL_Pele_2.

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
    new_index, _, original_dir = get_directory_new_index(pele_dir)

    if os.path.isdir(pele_dir):
        new_pele_dir = "{}_Pele_{}".format(original_dir, new_index)
        new_pele_dir = get_next_peledir(new_pele_dir)
        return os.path.join(os.path.dirname(pele_dir), new_pele_dir)
    else:
        return pele_dir


def get_latest_peledir(pele_dir):
    """
    Given a PELE directory it will return the name of the directory that
    looks newer. It employs the following criterion:

     - In case that NOL_Pele, NOL_Pele_1 and NOL_Pele_2 folders already
       exist in the working directory, it will choose NOL_Pele_2 since it
       has the highest suffix index.
     - In case no directory named NOL_Pele is found, the original name
       will be employed.

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
    new_index, old_index, original_dir = get_directory_new_index(pele_dir)

    # If the basic LIG_Pele already exists...
    if os.path.isdir(pele_dir):
        # If the potential "next" pele directory doesn't exist, the current pele_dir is the latest.
        latest_pele_dir = f"{original_dir}_Pele_{new_index}"
        if not os.path.isdir(latest_pele_dir):
            return pele_dir
        else:
            # Otherwise enumerate another one
            latest_pele_dir = get_latest_peledir(latest_pele_dir)
            return os.path.join(os.path.dirname(pele_dir), latest_pele_dir)
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
                        raise ValueError(
                            f"Check atom distance entrance {atom}. Should be like this: 'A:220:OD1'"
                        )
                    if (
                        line[21].strip() == chain
                        and line[22:26].strip() == resnum
                        and line[12:16].strip() == atomname
                    ):
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
    with open(pdb, "r") as f:
        waters = list(
            set(
                [
                    "{}:{}".format(line[21:22], line[22:26].strip())
                    for line in f
                    if line and "HOH" in line and (line.startswith("ATOM") or line.startswith("HETATM"))
                ]
            )
        )
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
            constraint = CONSTR_ATOM_POINT.format(
                spring_constant, chain, residue, atom_name
            )
        # Atom to atom constraint: 2.2-2.75-A:123:2-B:2:7 or 2.2-2.74-1985-1962
        elif len(constraint.split("-")) == 4:
            spring_constant, eq_distance, atom1_info, atom2_info = constraint.split("-")
            chain1, residue1, atom_name1 = retrieve_atom_info(atom1_info, pdb).split(
                ":"
            )
            chain2, residue2, atom_name2 = retrieve_atom_info(atom2_info, pdb).split(
                ":"
            )
            constraint = CONSTR_ATOM_ATOM.format(
                spring_constant,
                eq_distance,
                chain1,
                residue1,
                atom_name1,
                chain2,
                residue2,
                atom_name2,
            )
        final_constraints.append(constraint)
    return final_constraints


def retrieve_box(structure, residue_1, residue_2, weights=[0.5, 0.5]):
    # get center of interface (if PPI)
    coords1 = get_coords_from_residue(structure, residue_1)
    coords2 = get_coords_from_residue(structure, residue_2)

    box_center = np.average([coords1, coords2], axis=0, weights=weights)
    box_radius = (
        abs(np.linalg.norm(coords1 - coords2)) / 2 + 4
    )  # Sum 4 to give more space
    return list(box_center), box_radius


def get_coords_from_residue(structure, original_residue):
    parser = PDBParser()
    structure = parser.get_structure("protein", structure)
    chain, res_number, atom_name = original_residue.split(":")
    try:
        res_number = int(res_number)
    except ValueError:
        raise cs.WrongAtomStringFormat(
            f"The specified atom is wrong '{original_residue}'. \
Should be 'chain:resnumber:atomname'"
        )
    for residue in structure.get_residues():
        if residue.id[1] == res_number:
            for atom in residue.get_atoms():
                if atom.name == atom_name:
                    COI = np.array(list(atom.get_vector()))
                    return COI
    raise cs.WrongAtomSpecified(
        f"Atom {original_residue} could not be found in structure"
    )


def backup_logger(logger, message):
    if not logger:
        logger = logging.getLogger("logger")
        logger.setLevel(logging.INFO)
        logger.info(message)
    else:
        logger.info(message)


def find_nonstd_residue(pdb):
    with open(pdb, "r") as f:
        resnames = list(
            set(
                [
                    line[17:20]
                    for line in f
                    if line.startswith("ATOM")
                    and line[17:20] not in gv.default_supported_aminoacids
                ]
            )
        )
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
        raise ModuleNotFoundError(
            "Please install rdkit with the following command: conda install -c conda-forge rdkit"
        )


def get_suffix(filename, separator="_"):
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
    suffix = name.split("_")[-1]

    return suffix


def check_make_folder(output_folder):
    """
    Checks if output folders for plots exists and creates it, if not.

    Parameters
    ----------
    output_folder : str
        Name of the desired output folder
    """

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)


def check_remove_folder(*output_folders):
    """
    Removes the whole folder tree.

    Parameters
    ----------
    output_folders : Union[str, List[str]]
        Path(s) to folder to be removed
    """
    for folder in output_folders:
        if os.path.exists(folder):
            shutil.rmtree(folder, ignore_errors=True)


def get_atom_indices(ids, pdb, pdb_atom_name=None):
    """
    Extracts atom indices from a PDB file based on a tuple of chain IDs and residue names.

    Parameters
    -----------
    ids : List[Tuple(str, int)]
        List of atom IDs to be mapped
    pdb : str
        Path to PDB file.
    pdb_atom_name : str
        Set to desired PDB atom name, otherwise the whole residue will be extracted.

    Returns
    --------
    atom_indices : List[int]
        List of atom indices extracted from PDB
    """

    # Extract residue numbers and chain IDs from the list of tuples
    resnums = [element[1] for element in ids] if ids else []
    chains = [element[0] for element in ids] if ids else []

    # Extract relevant lines from PDB file
    with open(pdb, "r") as file:
        lines = file.readlines()
        lines = [
            line
            for line in lines
            if line.startswith("ATOM") or line.startswith("HETATM")
        ]

    atom_indices = list()

    for resnum, chain in zip(resnums, chains):
        for line in lines:
            # If chain and residue number match...
            if line[21:22].strip() == chain and line[22:26].strip() == str(resnum):
                atom_id = line[6:11].strip()

                # If the user specified PDB atom name, check that as well before appending
                if not pdb_atom_name:
                    atom_indices.append(atom_id)
                else:
                    if line[12:16].strip() == pdb_atom_name:
                        atom_indices.append(atom_id)

    # Explicitly convert to integers (mdtraj complains without dtype)
    atom_indices = list(
        np.array([int(atom_idx) for atom_idx in atom_indices], dtype=int)
    )
    return atom_indices


def retrieve_atom_names(pdb_file, residues):
    """
    Retrieves PDB atom names for a specific residue.

    Parameters
    -----------
    pdb_file : str
        Path to PDB file.
    residues : List[str]
        List of residue names for which the PDB atom names should be extracted.
    """
    from collections import defaultdict

    output = defaultdict(list)
    extracted_residues = list()

    with open(pdb_file, "r") as f:
        lines = f.readlines()

    # Retrieve HETATM lines only
    pdb_lines = [line for line in lines if line.startswith("HETATM") or line.startswith("ATOM")]

    for index, line in enumerate(pdb_lines):

        found_residue_name = line[17:20].strip()
        found_residue_number = line[22:26].strip()

        # Check if atom names for this residue are supposed to be extracted (or where already)
        if (
            found_residue_name in residues
            and found_residue_name not in extracted_residues
        ):
            atom_name = line[12:16]
            output[found_residue_name].append(atom_name)

        # Mark residue as extracted, if the next line has a different residue number
        try:
            next_residue_number = pdb_lines[index + 1][22:26].strip()
            if next_residue_number != found_residue_number:
                extracted_residues.append(found_residue_name)
        except IndexError:  # In case it's the end of PDB and index+1 doesn't exist
            continue

    return output


def get_residue_name(pdb_file, chain, residue_number):
    """
    Retrieves residue name from a PDB file based on residue number and chain.

    Parameters
    ------------
    pdb_file : str
        Path to PDB file.
    chain : str
        Chain ID.
    residue_number : str
        Number of the residue to check.

    Returns
    ---------
    resname : str
        Three letter name of the residue in lowercase, e.g. "cys".
    """
    if isinstance(residue_number, int):
        residue_number = str(residue_number)

    with open(pdb_file, "r") as file:
        lines = [line for line in file.readlines() if line.startswith("ATOM") or line.startswith("HETATM")]

    for line in lines:
        if line[21].strip() == chain and line[22:26].strip() == residue_number:
            residue_name = line[17:20].strip().lower()
            return residue_name
    else:
        raise ResidueNotFound(f"Could not find residue {residue_number} in chain {chain} in PDB file {pdb_file}.")


def get_residue_number(pdb_file, chain, residue_name):
    """
    Retrieves residue number from a PDB file based on chain ID and residue name.

    Parameters
    ----------
    pdb_file : str
        Path to PDB file.
    chain : str
        Chain ID of the residue.
    residue_name : str
        Residue name.

    Returns
    -------
        A string with residue number.
    """
    with open(pdb_file, "r") as file:
        pdb_lines = [line for line in file.readlines() if line.startswith("ATOM") or line.startswith("HETATM")]

    for line in pdb_lines:
        if line[21].strip() == chain and line[17:20].strip().lower() == residue_name.lower():
            residue_number = line[22:26].strip()
            return residue_number
    else:
        raise ResidueNotFound(f"Could not find residue {residue_name} in chain {chain} in PDB file {pdb_file}.")


# def get_pele_version():
#     """
#     Gets PELE version based on what's found under PELE variable.
#
#     Returns
#     -------
#     current_version : packaging.version
#         Version of PELE.
#
#     Raises
#     -------
#     PELENotFound if PELE variable is not set.
#     """
#     pele_path = os.path.join(pele_from_os, "bin/Pele_mpi")
#     output = subprocess.check_output(f"{pele_path} --version", shell=True)
#     current_version = version.parse(
#         re.findall(r"Version\: (\d+\.\d+\.\d+)\.", str(output))[0]
#     )
#
#     return current_version


def parse_atom_dist(atom_dist, pdb):
    """
    Parses atom distances defined by the user and checks their type.

    Parameters
    ----------
    atom_dist : str
        A single item of args.atom_dist set by the user.
    pdb : str
        Path to the PDB file (args.system).

    Returns
    -------
    atom : str
        Atom string in a format of chain:resnum:atom_name or residue string in the format of chain:resnum.
    atom_tag : str
        Tag necessary for the metrics JSON - "links" when using residue string and "atoms" when setting residue string.
    """
    if not str(atom_dist).isdigit() and len(atom_dist.split(":")) == 2:  # When atom_dist is a residue string
        atom = atom_dist
        atom_tag = "links"
    else:
        # When atom_dist is an atom number or an atom string
        atom = retrieve_atom_info(atom_dist, pdb)
        atom_tag = "atoms"
    return atom, atom_tag
