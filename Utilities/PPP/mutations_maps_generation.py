import sys

from MSM_PELE.PPP.program_own_classes import ZMATRIX

__author__ = 'jelisa'


def GenerateMap(zmat1, zmat2, name1, name2):
    # print name1, name2
    atoms_names1 = zmat1.AtomNames
    atoms_names2 = zmat2.AtomNames
    if len(atoms_names1) > len(atoms_names2):
        # print 'a'
        volatile_atoms_behaviour_ini = "disappear"
    elif len(atoms_names1) < len(atoms_names2):
        # print 'b'
        volatile_atoms_behaviour_ini = "appear"
    # else:
        # print "1. You're screwed, both aa. have the same number of atoms."
    common_atoms = []
    for at in atoms_names2:
        if at in atoms_names1:
            common_atoms.append(at)
        else:
            break
    # common_atoms = [at for at in atoms_names2 if at in atoms_names1]
    # print common_atoms
    if set(zmat1.GetBrotherAtoms(common_atoms[-1])) - set(common_atoms):
        last_common_parent = zmat1.GetParentName(common_atoms[-1])
        # print 'a', last_common_parent
        if last_common_parent not in common_atoms:
            print "u o something's terribly wrong!" \
                  "error in mutation from {} to {} ".format(name1, name2)
        res1_change = zmat1.GetAllChildrenAtoms(last_common_parent)
        res2_change = zmat2.GetAllChildrenAtoms(last_common_parent)
        for at in zmat1.GetDirectChildrenAtoms(last_common_parent):
            try:
                common_atoms.pop(common_atoms.index(at))
            except ValueError: pass
    else:
        res1_change = zmat1.GetAllChildrenAtoms(common_atoms[-1])
        res2_change = zmat2.GetAllChildrenAtoms(common_atoms[-1])
        # for at in zmat1.GetDirectChildrenAtoms(common_atoms):
        #     try:
        #         common_atoms.pop(common_atoms.index(at))
        #     except ValueError: pass
    if 'GLY' in [name1, name2]:
        if name1 == "GLY":
            changing_atoms = [['HA2', 'HA'], ['HA3', 'CB']]
            volatile_atoms_behaviour = "appear"
            volatile_atoms_behaviour_cont = "disappear"
            index = 7  # It's a magic number matching the number of atoms in the glycine
            volatile_atoms = atoms_names2[index:]
        else:
            changing_atoms = [['HA', 'HA2'], ['CB', 'HA3']]
            volatile_atoms_behaviour = "disappear"
            volatile_atoms_behaviour_cont = "appear"
            index = 7  # It's a magic number matching the number of atoms in the glycine
            volatile_atoms = atoms_names1[index:]
    else:
        common_changing_atoms = set(res1_change) & set(res2_change)
        changing_atoms = []
        for at in common_changing_atoms:
            changing_atoms.append([at, at])
            res1_change.pop(res1_change.index(at))
            res2_change.pop(res2_change.index(at))
        for at1, at2 in zip(res1_change, res2_change):
            changing_atoms.append([at1, at2])
            res1_change.pop(res1_change.index(at1))
            res2_change.pop(res2_change.index(at2))
        # print changing_atoms
        if len(atoms_names1) < len(atoms_names2):
            volatile_atoms = res2_change
            volatile_atoms_behaviour = "appear"
            volatile_atoms_behaviour_cont = "disappear"
        elif len(atoms_names1) > len(atoms_names2):
            volatile_atoms = res1_change
            volatile_atoms_behaviour = "disappear"
            volatile_atoms_behaviour_cont = "appear"
        else:
            volatile_atoms = []
            volatile_atoms_behaviour = volatile_atoms_behaviour_cont = ''
            # print "You're screwed, both aa. have the same number of atoms."
        # print 'bu'
    if "PRO" in [name1, name2]:
        if name1 == "PRO":
            changing_atoms.append(['CD', 'H'])
        else:
            changing_atoms.append(['H', 'CD'])
            try: volatile_atoms.pop(volatile_atoms.index('CD'))
            except ValueError:
                try:
                    common_atoms.pop(common_atoms.index('CD'))
                except ValueError:
                    print 'something is wrong with prolines...'
                    print name1, name2
                    print common_atoms
                    print changing_atoms
                    print volatile_atoms


    dir1 = "{}-{}".format(name1, name2)
    dir2 = "{}-{}".format(name2, name1)
    premap = [common_atoms, changing_atoms, volatile_atoms,
              {dir1: [volatile_atoms_behaviour, 0, 1],
               dir2: [volatile_atoms_behaviour_cont, 1, 0]}]

    dictio = {"{}-{}".format(name1, name2): premap}
    return dictio

# def main():
aminoacids = ['ALA', 'ARG', 'ASH', 'ASN', 'ASP', 'CYS', 'GLH', 'GLN', 'GLU',
              'GLY', 'HID', 'HIE', 'HIP', 'ILE', 'LEU', 'LYN', 'LYS', 'MET',
              'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

cyclics = ['HID', 'HIE', 'HIP', 'PHE', 'TRP', 'TYR', 'PRO']

# mutation_map = GenerateMap(ZMATRIX('ARG'), ZMATRIX('ASH'), 'ARG', 'ASH')
fileout = open(sys.argv[1], 'w')
fileout.write("{\n")
for x, aa1 in enumerate(aminoacids):
    for aa2 in aminoacids[x + 1:]:
        if aa1 == "PRO":
            print "mutation {} to {} do it by hand".format(aa1, aa2)
            continue
        elif aa1 in cyclics and aa2 in cyclics:
            print "mutation {} to {} do it by hand".format(aa1, aa2)
            continue
        elif aa1 == "GLY" and aa2 == "PRO":
            print "mutation {} to {} do it by hand".format(aa1, aa2)
            continue

        # print "working with :" , aa1, aa2
        zmap1 = ZMATRIX(aa1)
        zmap2 = ZMATRIX(aa2)
        mutation_map = GenerateMap(zmap1, zmap2, aa1, aa2)
        fileout.write("  '{}-{}':[\n".format(aa1, aa2))
        for x in mutation_map.values()[0]:
            fileout.write("     {},\n".format(x))
        # break
        fileout.write('  ],\n')
    # break

fileout.write('}')
fileout.close()

# if __name__ == '__init__':
#     main()
