import os

__author__ = 'jelisa'
DIR = os.path.dirname(__file__)

class ZMATRIX:
    """
    A class to store de data from the template file for PELE
    """

    AtomNames = None

    def __init__(self, aminoacid_name):
        # First the attributes are created as lists to be able to modify them when reading the data
        self.AtomNums = []
        self.Name = aminoacid_name
        self.AtomParents = []
        self.AtomNames = []
        self.BondLengths = []
        self.BondAngles = []
        self.Diherdrals = []
        self.VdWRadius = []
        self.ReadZmatrix(aminoacid_name)
        # The attributes are transformed to tuple so they can't be changed by mistake later on. It's really easy to
        # modify a list without noticing.
        self.AtomNums = tuple(self.AtomNums)
        self.AtomParents = tuple(self.AtomParents)
        self.AtomNames = tuple(self.AtomNames)
        self.BondLengths = tuple(self.BondLengths)
        self.BondAngles = tuple(self.BondAngles)
        self.Diherdrals = tuple(self.Diherdrals)
        self.VdWRadius = tuple(self.VdWRadius)

    def ReadZmatrix(self, aminoacid_name):
        """This function loads the zmatrix information for the desired aminoacid
        :rtype : a ZMATRIX object
        :param aminoacid_name: a string containing the 3 letters code for the desired aminoacid.
        """
        filename = "DeataLocal/Templates/OPLS2005/Protein/" + aminoacid_name.lower()
        try:
            filein = open(filename, 'r')
        except IOError:
            filename = "DaetaLocal/Templates/OPLS2005/HeteroAtoms/" + aminoacid_name.lower().strip() + 'z'
            try:
                filein = open(filename, 'r')
            except IOError:
                filename = os.path.join(DIR, "Data/Templates/OPLS2005/Protein/") + aminoacid_name.lower()
                try:
                    filein = open(filename, 'r')
                except IOError:
                    filename = os.path.join(DIR, "Data/Templates/OPLS2005/HeteroAtoms/") + aminoacid_name.lower() + 'z'
                    try:
                        filein = open(filename, 'r')
                    except IOError:
                        print("    * No such file or directory: {}".format(filename))
                        self.Name = None
                        return None
                    else:
                        self.Name = aminoacid_name.upper() + 'Z'
                else:
                    self.Name = aminoacid_name.upper()
            else:
                self.Name = aminoacid_name.upper() + 'Z'
        else:
            self.Name = aminoacid_name.upper()
        raw_data = filein.readlines()
        filein.close()
        add_lines = False
        add_vdw = False
        num_of_atoms = 0
        counter = 0
        for line in raw_data:
            if line.strip() == "":
                print("The template for the residue {} in {} has a blank line. " \
                      "Check it or not! Whatever you want.".format(self.Name, filename))
                continue
            if self.Name[-1] in ['z', 'b', 'e', 'a']:
                possible_name = line[:5]
            else:
                tmp_list = line.split()
                resname = tmp_list[0]
            if tmp_list[0] == self.Name:
                set_name = True
            else:
                possible_name = "".join(line[:5].strip().split())
                if possible_name == self.Name or possible_name == self.Name[:-1]:
                    set_name = True
                else:
                    set_name = False
            if set_name:
                num_of_atoms = int(tmp_list[1])
                counter = num_of_atoms
                add_lines = True
            elif add_lines:
                self.AtomNums.append(int(tmp_list[0]))
                self.AtomParents.append(int(tmp_list[1]))
                self.AtomNames.append(tmp_list[4].strip('_'))
                self.BondLengths.append(float(tmp_list[6]))
                self.BondAngles.append(float(tmp_list[7]))
                self.Diherdrals.append(float(tmp_list[8]))
                counter -= 1
                if counter == 0:
                    add_lines = False
            elif tmp_list[0] == "NBON":
                add_vdw = True
                counter = num_of_atoms
                if counter == 0:
                    print(" - There's a problem with the template of the residue {}.".format(aminoacid_name))
                    break
            elif add_vdw:
                self.VdWRadius.append((float(tmp_list[1]) / 2))
                counter -= 1
                if counter == 0:
                    break


    def GetParentName(self, atom2place):
        atom2place_position = self.AtomNames.index(atom2place)
        parent_num = self.AtomParents[atom2place_position]
        parent_position = self.AtomNums.index(parent_num)
        parent_name = self.AtomNames[parent_position]
        return parent_name

    def GetDirectChildrenAtoms(self, parent_atom_name):
        parent_atom_number = self.AtomNums[self.AtomNames.index(parent_atom_name)]
        if parent_atom_number in self.AtomParents:
            number_of_direct_children = self.AtomParents.count(parent_atom_number)
            children_indexes = [self.AtomParents.index(parent_atom_number)]
            for x in range(number_of_direct_children - 1):
                new_index = self.AtomParents[children_indexes[-1] + 1:].index(parent_atom_number)
                next_children = new_index + children_indexes[-1] + 1
                # The +1 is due to the fact that lists indexes start from 0
                children_indexes.append(next_children)
            children = [self.AtomNames[index] for index in children_indexes]
        else:
            children = []
        return children

    def GetAllChildrenAtoms(self, parent_atom_name, recursion_count=0):
        children = self.GetDirectChildrenAtoms(parent_atom_name)
        if recursion_count == 10000:
            raise RuntimeError("The function GetChildrenAtoms has arrived to the 10000 iterations!\n"
                               "             It's impossible that an aminoacid has so many atoms. Review it")
        for child in children:
            child_number = self.AtomNums[self.AtomNames.index(child)]
            if child_number in self.AtomParents:
                recursion_count += 1
                child_children = self.GetDirectChildrenAtoms(child)
                children.extend(child_children)
        return children

    def GetBrotherAtoms(self, atom_name):
        parent_name = self.GetParentName(atom_name)
        brothers = [brother for brother in self.GetDirectChildrenAtoms(parent_name)
                    if brother != atom_name and brother not in ['CA', 'C', 'N', 'O']]  # and brother[0] != 'H']
        return brothers

    def GetData2ComputeCoords(self, atom_name):
        at3_name = self.GetParentName(atom_name)
        at2_name = self.GetParentName(at3_name)
        at1_name = self.GetParentName(at2_name)
        atom_position_in_zmatrix = self.AtomNames.index(atom_name)
        r = self.BondLengths[atom_position_in_zmatrix]
        deta = self.BondAngles[atom_position_in_zmatrix]
        fi = self.Diherdrals[atom_position_in_zmatrix]
        return at1_name, at2_name, at3_name, r, deta, fi

    def ComputeDeltaFi(self):
        self.DeltaFi = [None for at in self.AtomNums]
        self.DeltaFiReference = [None for at in self.AtomNums]
        for index, atom in enumerate(self.AtomNames):
            if self.DeltaFi[index] is not None:
                continue
            elif atom in ["N", "H", "CA", "HA", "C", "O", "CB"]:
                self.DeltaFi[index] = 0
                continue
            else:
                brothers = self.GetBrotherAtoms(atom)
                if brothers:
                    if atom[0] != "H":
                        reference = atom
                    else:
                        for at in brothers:
                            if at[0] != "H":
                                reference = at
                                break
                        else:
                            reference = atom
                    self.DeltaFi[self.AtomNames.index(reference)] = 0
                    ref_fi = self.Diherdrals[self.AtomNames.index(reference)]
                    for at in brothers + [atom]:
                        if at == reference:
                            continue
                        else:
                            at_chi = self.Diherdrals[self.AtomNames.index(at)]
                            delta_fi = at_chi - ref_fi
                            self.DeltaFi[self.AtomNames.index(at)] = delta_fi
                            self.DeltaFiReference[self.AtomNames.index(at)] = reference
                else:
                    self.DeltaFi[index] = 0
        self.DeltaFi = tuple(self.DeltaFi)


class ROTAMERLIB:
    def __init__(self, aminoacid):
        self.Resname = ''
        self.NumberOfDihedrals = 0
        self.NumberOfPolarH = 0
        self.NumberOfDihedralAtoms = 0
        self.DihedralsDictio = {}
        self.PolarHydrogensDictio = {}
        self.ReadRotamerFile(aminoacid)

    def ReadRotamerFile(self, aminoacid):
        residue_rotamer_library_path = os.path.join(DIR, "Data/RotamerLibs/") + aminoacid.upper() + '.side'
        filein = open(residue_rotamer_library_path, 'r')
        file_text = filein.readlines()
        filein.close()
        rotamer_info = file_text[0].split()
        self.Resname = rotamer_info[1]
        self.NumberOfDihedralAtoms = int(rotamer_info[2])
        self.NumberOfDihedrals = int(rotamer_info[3])
        multiplier = float(rotamer_info[4])
        self.NumberOfPolarH = int(rotamer_info[5])
        ordered_keys = [element[:-1].strip('_') for element in file_text[1:self.NumberOfDihedralAtoms + 1]]
        self.DihedralsDictio = {key: tuple([float(x.strip().split()[y]) * multiplier
                                            for x in file_text[self.NumberOfDihedralAtoms + self.NumberOfPolarH + 1:]])
                                for key, y in zip(ordered_keys, range(len(ordered_keys)))}
