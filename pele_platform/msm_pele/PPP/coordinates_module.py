import sys

from numpy import asarray, vstack, cos, cross, dot, pi, sin, transpose, linalg, arctan2


__author__ = 'jelisa'


fixed_atoms = ["N", "H", "CA", "HA", "C", "O", "CB"]


def Grades2Radians(angle):
    """This is a small function to convert from degrees to radians a given angle.
    :rtype : a float with the angle in radians
    :param angle: an angle in degrees
    :return: the angle in radians
    """
    return (angle / 180.) * pi


def Radians2Degrees(angle):
    """
    :rtype : float
    :param angle: float with the angle in radians
    :return: the angle in degrees
    """
    return (angle * 180) / pi


def ComputeAngleDifference(angle1, angle2):
    """
    This function computes the difference between two angles.
    :rtype : float
    :param angle1: float with the first angle
    :param angle2: float with the second angle
    :return: 
    """
    if angle1 < 0:
        inicial = angle1 + 360
    else:
        inicial = angle1
    if angle2 < 0:
        final = angle2 + 360
    else:
        final = angle2
    diff = inicial - final
    return diff


# noinspection PyTypeChecker
def ComputeDihedral(at1, at2, at3, at4, v=False):
    """
    This function implements the calculation of the Dihedral using the arctan as explained
    in the link: http://math.stackexchange.com/questions/47059/how-do-i-calculate-a-dihedral-angle-given-cartesian-coordinates
    date 14/05/2015. It's the same algorithm as the one implemented in PELE-1.5++ at date 15/05/2015
    :rtype : float
    :param at1: coordinates of the first atom
    :param at2: coordinates of the second atom
    :param at3: coordinates of the third atom
    :param at4: coordinates of the fourth atom
    :return: the dihedral angle in degrees between the four atoms.
    """
    b1 = at1 - at2
    b2 = at2 - at3
    b3 = at3 - at4
    c_b1b2 = cross(b1, b2)
    c_b2b3 = cross(b2, b3)
    u_b2 = b2 / linalg.norm(b2)
    # print u_b2
    n1 = c_b1b2 / linalg.norm(c_b1b2)
    n2 = c_b2b3 / linalg.norm(c_b2b3)
    # print "compute dihedrals, a", n1.shape, n2.shape
    m1 = cross(n1, u_b2)
    x = dot(n1, n2)
    y = dot(m1, n2)
    z = dot(u_b2, n2)
    # print z
    if z > 0.0001:
        print("Z should be almost 0, something in wrong!")
        sys.exit()
    dihedral = arctan2(y, x)
    if v:
        print("dihedral in rad:", dihedral)
    return Radians2Degrees(dihedral)


# noinspection PyTypeChecker
def ComputeCartesianCoordinates(at1, at2, at3, r, deta, fi):
    """
    This function implements the SN-NeRF algorithm described in Parson et al., 2005  (Ref. 14)
    :rtype : array with the coordinates
    :param at1: An AtomGroup from prody containing the atom 1 (atom A in the article)
    :param at2: An AtomGroup from prody containing the atom 2 (atom B in the article)
    :param at3: An AtomGroup from prody containing the atom 3 (atom C in the article)
    :param r: The distance between the atom to place and at3.
    :param deta: The angle defined by at2 and the atom to place with the at3 as vertex
    :param fi: The dihedral angle defined by at1, at2, at3 and the atom to place
    """

    rad_deta = Grades2Radians(deta)
    rad_fi = Grades2Radians(fi)
    # print 'compute coords', fi, rad_fi
    bc = at3.getCoords()[0] - at2.getCoords()[0]
    ab = at2.getCoords()[0] - at1.getCoords()[0]
    u_bc = bc / linalg.norm(bc)
    n = cross(ab, u_bc) / linalg.norm(cross(ab, u_bc))
    m = transpose(vstack([u_bc, cross(n, u_bc), n]))
    d2 = asarray([-r * cos(rad_deta), r * cos(rad_fi) * sin(rad_deta), r * sin(rad_fi) * sin(rad_deta)])
    coords = dot(m, d2) + at3.getCoords()[0]
    return asarray(coords)


def GenerateCoordinatesFromZmatrix(initial_residue, atoms2add, zmatrix):
    atoms2add_coords = []
    for atom in atoms2add:
        at1_name, at2_name, at3_name, r, deta, fi = zmatrix.GetData2ComputeCoords(atom)
        # print atom, at1_name, at2_name, at3_name
        # print zmatrix.Name, zmatrix.AtomNames
        at3 = initial_residue.select('name {}'.format(at3_name))
        at2 = initial_residue.select('name {}'.format(at2_name))
        at1 = initial_residue.select('name {}'.format(at1_name))
        # print initial_residue.getNames()
        # print at1, at2, at3
        atoms2add_coords.append(ComputeCartesianCoordinates(at1, at2, at3, r, deta, fi))
    return asarray(atoms2add_coords)


def ChangeResidueCoordinates(initial_residue, zmatrix, rotamer_library, dihedral_index):
    """
    This function changes the coordinates of a residue according to the
    information in the zmatrix and the rotamer library used in PELE.

    :type rotamer_library: ROTAMERLIB
    :type zmatrix: ZMATRIX
    :type initial_residue: AtomGroup of prody
    :type dihedral_index: int
    :rtype : AtomGroup of prody
    :param initial_residue: atom group of prody with the initial residue to modify
    :param zmatrix: ZMATRIX object with the information for the residue
    :param rotamer_library: ROTAMERLIB object containing the information of the rotamer library for the residue
    :param dihedral_index: the index that indicates which rotamer should be used
    :return: AtomGroup of prody with the coordinates changed
    """
    residue = initial_residue.copy()
    dihedral_atoms_brothers = []
    for atom_name in rotamer_library.DihedralsDictio.keys():
        dihedral_atoms_brothers.extend(zmatrix.GetBrotherAtoms(atom_name))
    atoms_yet_to_place = []
    for atom_name in zmatrix.AtomNames:
        if atom_name in fixed_atoms:
            continue
        elif atom_name[0] == "H":
            atoms_yet_to_place.append(atom_name)
            continue
        else:
            at1_name, at2_name, at3_name, r, deta, zmatrix_fi = zmatrix.GetData2ComputeCoords(atom_name)
            at3 = residue.select("name {}".format(at3_name))
            at2 = residue.select("name {}".format(at2_name))
            at1 = residue.select("name {}".format(at1_name))
            atom = residue.select("name {}".format(atom_name))
            delta_fi = zmatrix.DeltaFi[zmatrix.AtomNames.index(atom_name)]
            if atom_name in rotamer_library.DihedralsDictio.keys():
                initial_fi = rotamer_library.DihedralsDictio[atom_name][dihedral_index]
            elif atom_name in dihedral_atoms_brothers:
                dihedral_brother = [at for at in zmatrix.GetBrotherAtoms(atom_name)
                                    if at in rotamer_library.DihedralsDictio.keys()]
                try:
                    dihedral_brother[0]
                except IndexError:
                    print("Something went wrong when checking brothers while moving the residue to solve clashes.")
                    sys.exit()
                else:
                    initial_fi = ComputeDihedral(at1.getCoords()[0], at2.getCoords()[0], at3.getCoords()[0],
                                                 residue.select("name {}".format(dihedral_brother[0])).getCoords()[0])
            elif zmatrix.DeltaFi[zmatrix.AtomNames.index(atom_name)] != 0 and zmatrix.GetBrotherAtoms(atom_name):
                brother = zmatrix.GetBrotherAtoms(atom_name)[0]
                initial_fi = ComputeDihedral(at1.getCoords()[0], at2.getCoords()[0], at3.getCoords()[0],
                                             residue.select("name {}".format(brother)).getCoords()[0])
            else:
                initial_fi = zmatrix_fi

            fi = initial_fi + delta_fi
            coords = ComputeCartesianCoordinates(at1, at2, at3, r, deta, fi)
            atom.setCoords(coords)
    for atom_name in atoms_yet_to_place:
        at1_name, at2_name, at3_name, r, deta, zmatrix_fi = zmatrix.GetData2ComputeCoords(atom_name)
        at3 = residue.select('name {}'.format(at3_name))
        at2 = residue.select('name {}'.format(at2_name))
        at1 = residue.select('name {}'.format(at1_name))
        atom = residue.select('name {}'.format(atom_name))
        delta_fi = zmatrix.DeltaFi[zmatrix.AtomNames.index(atom_name)]
        if atom_name in dihedral_atoms_brothers:
            dihedral_brother = [at for at in zmatrix.GetBrotherAtoms(atom_name)
                                if at in rotamer_library.DihedralsDictio.keys()]
            try:
                dihedral_brother[0]
            except IndexError:
                print("Something went wrong when checking brothers while moving the residue to solve clashes.")
                sys.exit()
            else:
                initial_fi = ComputeDihedral(at1.getCoords()[0], at2.getCoords()[0], at3.getCoords()[0],
                                             residue.select('name {}'.format(dihedral_brother[0])).getCoords()[0])
        elif zmatrix.DeltaFi[zmatrix.AtomNames.index(atom_name)] != 0 and zmatrix.GetBrotherAtoms(atom_name):
            brother = zmatrix.GetBrotherAtoms(atom_name)[0]
            initial_fi = ComputeDihedral(at1.getCoords()[0], at2.getCoords()[0], at3.getCoords()[0],
                                         residue.select("name {}".format(brother)).getCoords()[0])
        else:
            initial_fi = zmatrix_fi
        fi = initial_fi + delta_fi
        coords = ComputeCartesianCoordinates(at1, at2, at3, r, deta, fi)
        atom.setCoords(coords)
    return residue


def ModifyCoords(initial_residue, zmatrix, atom_name):
    zmatrix.ComputeDeltaFi()
    at1 = initial_residue.select("name N")
    at2 = initial_residue.select("name C")
    at3 = initial_residue.select("name CA")
    atom = initial_residue.select("name {}".format(atom_name))
    r = zmatrix.BondLengths[zmatrix.AtomNames.index(atom_name)]
    deta = zmatrix.BondAngles[zmatrix.AtomNames.index(atom_name)]
    fi = ComputeDihedral(at1.getCoords()[0], at2.getCoords()[0], at3.getCoords()[0], atom.getCoords()[0])
    coords = ComputeCartesianCoordinates(at1, at2, at3, r, deta, fi)
    return coords


def ModifyCoordinatesPRO(residue, zmatrix, atom_name, reference_atom):
    parent_name = zmatrix.GetParentName(atom_name)
    parent_atom = residue.select("name {}".format(parent_name)).copy()
    atom = residue.select("name {}".format(reference_atom)).copy()
    parent_atom_coords = parent_atom.getCoords()
    vector = atom.getCoords() - parent_atom_coords
    unitary_vector = vector / linalg.norm(vector)
    r = zmatrix.BondLengths[zmatrix.AtomNames.index(atom_name)]
    atom_coords = parent_atom_coords + unitary_vector * r
    return atom_coords
