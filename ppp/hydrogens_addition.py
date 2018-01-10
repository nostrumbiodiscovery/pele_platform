from prody import AtomGroup

from mutational_module import AddAtoms, ModifyExistingAtoms
from program_own_classes import ZMATRIX
from coordinates_module import ComputeCartesianCoordinates, ComputeDihedral
from global_processes import FindInitialAndFinalResidues

__author__ = 'jelisa'


def DefineNewAtom(atom_name, element, coordinates, resname, resnum, chain_id):
    """
    This function creates a new AtomGroup instance containing one atom.
    """
    new_atom = AtomGroup()
    new_atom.setNames([atom_name])
    new_atom.setElements([element])
    new_atom.setCoords([coordinates])
    new_atom.setResnames(resname)
    new_atom.setResnums(resnum)
    new_atom.setChids(chain_id)
    new_atom.setAltlocs([''])
    new_atom.setBetas([0])
    new_atom.setIcodes([''])
    new_atom.setOccupancies([1])
    new_atom.setSegnames([''])
    new_atom.setSerials([0])
    # new_atom.setAnisous([[0.0, 0.0, 0.0]])
    return new_atom


def PlaceSpecialAtoms(old_residue, atom_name, structure, resnum, zmatrix, verbose=False):
    current_residue = structure.select("resnum {}".format(resnum)).copy()
    atom_position_in_zmatrix = zmatrix.AtomNames.index(atom_name)
    r = zmatrix.BondLengths[atom_position_in_zmatrix]
    deta = zmatrix.BondAngles[atom_position_in_zmatrix]
    element = "H"
    if atom_name == "H":
        previous_residue = structure.select("resnum {}".format(resnum - 1))
        at3 = current_residue.select("name N")
        at2 = previous_residue.select("name C")
        at1 = previous_residue.select("name O")
        fi = 180
    elif atom_name in ["HA2", "HA"]:
        # print '2'
        previous_residue = structure.select("resnum {}".format(resnum - 1))
        if previous_residue is None:
            at4 = current_residue.select("name N")
            at3 = current_residue.select("name CA")
            at2 = current_residue.select("name C")
            at1 = current_residue.select("name O")
            fi = -120 + ComputeDihedral(at1.getCoords()[0], at2.getCoords()[0],
                                        at3.getCoords()[0], at4.getCoords()[0])
        else:
            at4 = current_residue.select("name C")
            at3 = current_residue.select("name CA")
            at2 = current_residue.select("name N")
            at1 = previous_residue.select("name C")
            fi = 120 + ComputeDihedral(at1.getCoords()[0], at2.getCoords()[0],
                                       at3.getCoords()[0], at4.getCoords()[0])
    elif atom_name in ["HA3", "CB"]:
        previous_residue = structure.select("resnum {}".format(resnum - 1))
        if previous_residue is None:
            at4 = current_residue.select("name N")
            at3 = current_residue.select("name CA")
            at2 = current_residue.select("name C")
            at1 = current_residue.select("name O")
            fi = -240 + ComputeDihedral(at1.getCoords()[0], at2.getCoords()[0],
                                        at3.getCoords()[0], at4.getCoords()[0])
        else:
            at4 = current_residue.select("name C")
            at3 = current_residue.select("name CA")
            at2 = current_residue.select("name N")
            at1 = previous_residue.select("name C")
            fi = 240 + ComputeDihedral(at1.getCoords()[0], at2.getCoords()[0],
                                       at3.getCoords()[0], at4.getCoords()[0])
    elif atom_name == "OXT":
        at4 = current_residue.select("name O")
        at3 = current_residue.select("name C")
        at2 = current_residue.select("name CA")
        at1 = current_residue.select("name N")
        fi = 180 + ComputeDihedral(at1.getCoords()[0], at2.getCoords()[0],
                                   at3.getCoords()[0], at4.getCoords()[0])
        element = "O"
    elif atom_name in ["H1", "H2", "H3"]:
        at3 = current_residue.select("name N")
        at2 = current_residue.select("name CA")
        if current_residue.getResnames()[0] == "GLY":
            at1 = current_residue.select("name HA2")
        else:
            at1 = current_residue.select("name CB")
        # This magic numbers come from placing the hydrogens in maestro
        # and computing the angle and the dihedrals.
        fi = 60 + 120 * (int(atom_name[1]) - 1)  # dihedral between C-Calpha-N-hydrogens in maestro
        deta = 109.53580  # Angle between Calpha-N-Hydrogens in maestro.
    else:
        print "How did you get here??"
        return old_residue
    coords = ComputeCartesianCoordinates(at1, at2, at3, r, deta, fi)
    new_atom = DefineNewAtom(atom_name, element, coords, current_residue.getResnames()[:1],
                             current_residue.getResnums()[:1], current_residue.getChids()[:1])
    # print len(old_residue.getAnisous()), len(new_atom.getAnisous())
    residue = old_residue + new_atom
    return residue


def PlaceHydrogen(old_residue, atom_name, zmatrix):
    at1_name, at2_name, at3_name, r, deta, zmatrix_fi = zmatrix.GetData2ComputeCoords(atom_name)
    # print 'a', old_residue.getNames(), old_residue.getResnames()
    # print 'b', zmatrix.AtomNames
    at1 = old_residue.select("name {}".format(at1_name))
    at2 = old_residue.select("name {}".format(at2_name))
    at3 = old_residue.select("name {}".format(at3_name))
    delta_fi = zmatrix.DeltaFi[zmatrix.AtomNames.index(atom_name)]
    reference_atom_name = zmatrix.DeltaFiReference[zmatrix.AtomNames.index(atom_name)]
    reference_atom = old_residue.select("name {}".format(reference_atom_name))
    # print 'd', zmatrix.DeltaFi
    # print 'c', at1, at2, at3, reference_atom_name, atom_name
    fi = ComputeDihedral(at1.getCoords()[0], at2.getCoords()[0], at3.getCoords()[0], reference_atom.getCoords()[0]) + delta_fi
    coords = ComputeCartesianCoordinates(at1, at2, at3, r, deta, fi)
    new_atom = DefineNewAtom(atom_name, "H", coords, old_residue.getResnames()[:1],
                             old_residue.getResnums()[:1], old_residue.getChids()[:1])
    # print old_residue.getAnisous(), new_atom.getAnisou()
    new_residue = old_residue + new_atom
    return new_residue


def FixStructure(initial_structure, residues2fix, debug=False):
    final_structure = None
    for chain in initial_structure.iterChains():
        chain_id = chain.getChid()
        if chain_id == " ":
            current_structure = initial_structure.copy()
        else:
            current_structure = initial_structure.select("chain `{}`".format(chain_id)).copy()
        new_chain = None
        initial_residue, final_residue = FindInitialAndFinalResidues(current_structure)
        for residue in current_structure.iterResidues():
            resId = " ".join([residue.getResname(), str(residue.getResnum()), residue.getChid()])
            old_res = current_structure.select("resname `{}` and resnum `{}`".format(residue.getResname(),
                                                                                     residue.getResnum())).copy()
            if resId not in residues2fix.keys():
                new_residue = old_res
            else:
                print " Adding to the residue '{}' the following atoms\n  {}".format(resId, residues2fix[resId])
                if residue.getResnum() == initial_residue:
                    zmatrix = ZMATRIX(residue.getResname() + 'B')
                elif residue.getResnum() == final_residue:
                    zmatrix = ZMATRIX(residue.getResname() + 'E')
                else:
                    try:
                        residues2fix[resId]
                    except KeyError:
                        zmatrix = ZMATRIX(residue.getResname())
                    else:
                        if "H2" in residues2fix[resId]:
                            zmatrix = ZMATRIX(residue.getResname() + 'B')
                        elif "OXT" in residues2fix[resId]:
                            zmatrix = ZMATRIX(residue.getResname() + 'E')
                        else:
                            zmatrix = ZMATRIX(residue.getResname())
                if zmatrix.Name is None:
                    print "The residue {} {} doesn't have a template, so it cannot be fixed.\n" \
                          "PELE won't work! Check it!".format(residue.getResname(), residue.getResnum())
                    continue
                if debug:
                    print "ZMATRIX used:", zmatrix.Name
                if residue.getResnum() == initial_residue:
                    atoms2add = {"H1", "H2", "H3"}.union(residues2fix[resId]).difference(set(residue.getNames()))
                elif residue.getResnum() == final_residue:
                    atoms2add = {"OXT"}.union(residues2fix[resId]).difference(set(residue.getNames()))
                else:
                    atoms2add = residues2fix[resId]
                if debug:
                    print 'atoms to add:', atoms2add
                    print 'residues2fix:', residues2fix[resId]
                residue_info = {"fin_resname": residue.getResname(),
                                "resnum": residue.getResnum(), "chain": residue.getChid()}
                if residue.getResname() == "HOH":
                    print zmatrix.Name
                    new_residue = current_structure.select('resnum {}'.format(residue.getResnum())).copy()
                else:
                    try:
                        zmatrix.ComputeDeltaFi()
                    except ValueError:
                        print "Something went wrong! Residue {} {}".format(residue.getResname(), residue.getResnum())
                    for name in zmatrix.AtomNames:
                        if name in atoms2add:
                            if name == "OXT" and "HXT" in residue.getNames():
                                atomnames_of_2_letters = ["FE"]
                                new_residue = ModifyExistingAtoms(residue, [["HXT", "OXT"]], atomnames_of_2_letters, 0, 1, zmatrix)
                                old_res = new_residue
                            elif name in ["H", "HA2", "HA", "HA3", "CB", "OXT", "H1", "H2", "H3"]:
                                if debug:
                                    print '000'
                                    new_residue = PlaceSpecialAtoms(old_res, name, current_structure, residue.getResnum(),
                                                                    zmatrix, True)
                                else:
                                    new_residue = PlaceSpecialAtoms(old_res, name, current_structure, residue.getResnum(),
                                                                    zmatrix)
                                old_res = new_residue
                            elif name[0] == "H" and zmatrix.DeltaFi[zmatrix.AtomNames.index(name)] != 0:
                                new_residue = PlaceHydrogen(old_res, name, zmatrix)
                                old_res = new_residue
                            else:
                                new_residue = AddAtoms(old_res, [name], ["FE"], residue_info, zmatrix, False)
                                old_res = new_residue
            if new_residue.getAnisous() is not None and not len(new_residue.getAnisous()) == new_residue.numAtoms():
                """If this values don't match this will give problems later when copying the structure.
                The magic number 6 used to define the new_anisous is due to the fact that the anisotropic
                arrays have six positions."""
                new_anisous = residue.getAnisous()
                new_anisous.resize(new_residue.numAtoms(), 6)
                new_residue.setAnisous(new_anisous)
                if debug:
                    print "checking the concordance between the anisou and the number of atoms"
                    print "len anisous: {}, number of atoms: {}".format(len(new_residue.getAnisous()),
                                                                        new_residue.numAtoms())
            if new_chain is None:
                new_chain = new_residue
            else:
                new_chain = new_chain + new_residue
        if final_structure is None:
            final_structure = new_chain
        else:
            final_structure = final_structure + new_chain
    final_structure.setTitle("Structure with Hydrogens")
    return final_structure
