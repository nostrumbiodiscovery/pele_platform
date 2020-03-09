import os
import sys
import re
import logging
import warnings

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


def preproces_lines(lines):
    for i, line in enumerate(lines):
        line = re.sub(' +', ' ', line)
        line = line.strip('\n').strip().split()
        lines[i] = line
    return lines


def is_repited(pele_dir):

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
                new_pele_dir = is_repited(new_pele_dir)
                return new_pele_dir
    else:
                return pele_dir


def is_last(pele_dir):

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
                    chain, resnum, atomname = atom.split(":")
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
        sys.exit("Check the atoms given to calculate the distance metric")


def find_coords(pdb, resnum, chain, atom="OW"):
    with open(pdb, "r") as f:
        for line in f:
            if line:
                if line.startswith("HETATM") and line[22:26].strip() == resnum and line[21:22] == chain:
                    return [float(coord) for coord in line[30:54].split()]


def find_centroid(points):
    x = [cx for cx, cy, cz in points]
    y = [cy for cx, cy, cz in points]
    z = [cz for cx, cy, cz in points]
    n_points = len(points)
    centroid = (sum(x) / n_points, sum(y) / n_points, sum(z) / n_points)
    return centroid

def retrieve_all_waters(pdb):
    with open(pdb, 'r') as f:
        return list(set(["{}:{}".format(line[21:22], line[23:26].strip()) for line in f if line and "HOH" in line]))

def retrieve_constraints_for_pele(constraints, pdb):
    CONSTR_ATOM_POINT = '{{ "type": "constrainAtomToPosition", "springConstant": {}, "equilibriumDistance": 0.0, "constrainThisAtom": "{}:{}:{}" }},'
    CONSTR_ATOM_ATOM = '{{"type": "constrainAtomsDistance", "springConstant": {}, "equilibriumDistance": {}, "constrainThisAtom":  "{}:{}:{}", "toThisOtherAtom": "{}:{}:{}"}},'
    final_constraints = []
    for constraint in constraints:
        #Atom to point constraint: 2.2-A:123:2 or 2.2-1986
        if len(constraint.split("-")) == 2:
            spring_constant, atom_info = constraint.split("-")
            chain, residue, atom_name = retrieve_atom_info(atom_info, pdb).split(":")
            constraint = CONSTR_ATOM_POINT.format(spring_constant, chain, residue, atom_name)
        #Atom to atom constraint: 2.2-2.75-A:123:2-B:2:7 or 2.2-2.74-1985-1962
        elif len(constraint.split("-")) == 4:
            spring_constant, eq_distance, atom1_info, atom2_info = constraint.split("-")
            chain1, residue1, atom_name1 = retrieve_atom_info(atom1_info, pdb).split(":")
            chain2, residue2, atom_name2 = retrieve_atom_info(atom2_info, pdb).split(":")
            constraint =  CONSTR_ATOM_ATOM.format(spring_constant, eq_distance, chain1, residue1, atom_name1, chain2, residue2, atom_name2)
        final_constraints.append(constraint)
    return final_constraints
        
