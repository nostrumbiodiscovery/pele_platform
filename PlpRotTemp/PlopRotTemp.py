about = """
$Revision: 1.13 $
For ligands:
Reads in a maestro mae file and makes a "rotamer enabled" template and the rotamer library to accompany it.  This consists of finding the backbone core that results in the least number of child bonds rotated with any rotatable bond rotation.  Reads in the rotatable bonds from a macromodel atomtyping (should be easily modifyiable to read them from stdin).  Hetgrp_ffgen is used for atomtyping and to determine the geometry in the template fromt the mae.  The mae must only have the residue to atomtype in it. 

Builds rotamer libraries for arbitrary ligand molecules by building and combining rotamer libraries.  There are two general algorithms that are implemented.  This first is using macromodel (or some other external tool) to sample the degrees of freedom and converting the resulting ensemble into a rotamer library.  The second is identifying the rotatable bonds, again using macromodel, and assigning rotamer libraries to these bonds.  For most bonds a simple freely rotatable library ( 0,10,20,30...360), but in the case of closed rings special libraries are built using macromodel sampling.  These component rotamer libraries are then arranged into groups for use in PLOP.  Each group consists of a single tree rooted at the central core.  This core can either be used chosen, or will be chosen based on an algorithm that minimizes the number of bond lengths from the farthest leeaf to the trunk.  Any built rotamer libraries are placed in the current directory and a file called <resname>.rot.assign is also written.  This tells PLOP how to assemble the full combinatoric library that will be used in sidehchain prediction/monte carlo.     

For unnatural amino acids:
Requires: 1) a maestro mae file of an unnatural amino acid with no NMA or ACE termini; the N-H and C=0 are to be left as they are found in a peptide
          2) a template file (specified by -t=<FILENAME>) created by hetgrp-ffgen using a maestro mae file of an unnatural amino acid with NMA or ACE termini present
Option Settings Required:  -unnat=yes -t=<FILENAME> [other options] <RESIDUE_WITHOUT_CAPPING_GROUPS>.mae
Outputs: 1) a re-ordered template file (old one is backed up in FILENMAE.hetgrp_ffgen)
         2) a PLOP nonstandard residue specification for pasting into a PLOP control file, both to stdout and to <maefile>_torsions.txt


Options
   -core=<an1>      Give one atom of the core section
   -f=<an1>,<an2>   Fix this bond number in the mae file.  
                    CANNOT USE MACORMODEL SAMPLING WITH THIS OPTION 
   -tor=<an1>,<an2>      Define torsions instead of using macromodel to do it.  
                         CANNOT USE MACORMODEL SAMPLING WITH THIS OPTION 
                         YOU MUST ALSO SPECIFY WHETHER THERE ARE RINGS OR NOT (-r) 
   -btor=<an1>,<an2>     Define torsions for the backbone  
   -m=<yes/no>           Wheter  to split into multiple libraries
   -min=<yes/no>         Whether to minimize
   -r=<yes/no/only>      Whether to allow flexible rings
   -t=<FILENAME>         Use this original template instead of creating one
   -o=<FILENAME>         Write output to this filename (default is the residue name)
   -g=<GRIDRES>          Input the grid resolution for plop rotamer library
   -goh=<GRIDRES>        Input grid resolution for plop rotamer library (OHs)
   -OPLS=<2001,2005>     OPLS version to use
   -mae_charges=<yes/no> Use mae charges for template file
   -mtor=<number>        Gives the maximum number of torsions allowed in each
                         group.  Will freeze bonds to extend the core if 
                         necessary.
   -clean=<yes/no>       Whether to clean up all the intermediate files
   -d                    Debuging mode, assume all external programs have been
                         already run and the intermediate files were not cleaned

Options for Macromodeling Sampling 
   seting any of these options will do the sampling in macromodel to build 
   a sidechain library.  Otherwise a series of libraries will be indentified 
   to be combined together in plop later.  

   -c=<FILENAME>    Use the following pdb/mae file as an exhaustive list of
                    all confomers.  This must be one file with multiple
                    models/entries.  Still uses macromodel to find rotatable
                    bonds. 
   -n=<number>      Maximum Number of Entries in Rotamer File
   -a=<CGEN,MCMM...> Type of Search to Run for Sidechains, if none is specified 
                    it will use PLOP sampling
   -ba=<CGEN,MCMM...>Type of Search to Run for Backbones, if none is specified it 
                    will not be sampled and a backbone library and no backbone
                    library will be created
   -s=<number>      Maximum Conformers to Sample
   -e=<number>      Energy Cutoff (in kJ/mole)
   -md=<number>     Maximum distance between atoms in equivalent struct(A)
Options for Unnatural Amino Acids
   -chain=<chain>   For output of PLOP-style nonstandard residue specification,
                    set the chain name equal to what it would be in the whole
                    macromolecule
   -res=<num>       For output of PLOP-style nonstandard residue specification,
                    set the residue number equal to what it would be in the whole
                    macromolecule

Mae file should be properly atomtyped

Most common problem:  As part of this procedure the pdb atom names are often renamed when they are not unique. Alsothis procedure works best if the ligand is minimized.  For these two reasons an atomtyped, minimzed version of the input ligand is written to (input).PlopRotTemp.pdb.  If at all possible, use the ligand structure and atom names from this file in any subsequent plop runs.    

examples:
Build a rotamer library for the following ligand at a grid resolution of 20 deg using PLOP/PRIME to combine libraries
$SCHRODINGER/utilities/python PlopRotTemp.py 3ert_lig.mae -g=20

Build a rotamer library for the following ligand at using CGEN sampling in macromodel.
$SCHRODINGER/utilities/python PlopRotTemp.py 3ert_lig.mae -a=CGEN

Build a rotamer library for the following ligand at a grid resolution of 20 deg using macromodel to sample any rings and combining this with freely rotatable libraries in PLOP to create combined libraries for the ligand.
$SCHRODINGER/utilities/python PlopRotTemp.py 1rth_lig.mae -r=yes

Make libraries for rotatable bonds in ligand.mae up to a maximum of 4 rotatable bonds in each library
All additional bonds are defined as backbone and are sampled with CGEN to produce a backbone libary
$SCHRODINGER/utilities/python PlopRotTemp.py ligand.mae -mtor=4 -ba=CGEN

For a given ligand named LIG the following files will be created:
lig                - Template file for use in PLOP, its zmatrix matches the libraries created
LIG.rot.assign     - Summary of all libraries build or used for this ligand read into plop with the command
                     "rot assign all"
LIG???.side        - (OPTIONAL) Component sidechains libraries created if there are closed rings or CGEN sampling is used
LIG__B.back        - (OPTIONAL) Backbone sidechain library

---------------------------------------------------------------------------------------------------------



All jobs run on the localhost

"""
################################################################################
# Globals/Constants 
################################################################################
usage = """
usage: "$SCHRODINGER/utilities/python PlopRotTemp.py [file.mae]"
 """

################################################################################
# Packages 
################################################################################
import sys
import re
import math
import array
import os
import warnings
import shutil
import subprocess
from rdkit import Chem
from rdkit.Chem import TorsionFingerprints
from schrodinger import structure
# from schrodinger import structureutil
hetgrp_ffgen = os.environ['SCHRODINGER'] + "/utilities/hetgrp_ffgen"



################################################################################
# Definitions 
################################################################################

ERROR_ATOMTYPES = 'ATOM NAMES REPITED IN MAE FILE'
ERROR_ROTAMER_LIB = 'ERROR: LACK OF NON BONDED OR BOND PRAMETERS IN ROTAMER LIGAND'
DEFAULT_RADIUS_VDW = '0.5000'
DEFAULT_SPRING_K = '268.0'
DEFAULT_EQ_DIST = '1.529'
FILE_DIR_PATH = os.path.dirname(__file__)
PARAM_PATH = os.path.join(FILE_DIR_PATH, 'param/sgbnp.param') ##impact??
SIMILARITY_PATH = os.path.join(FILE_DIR_PATH, 'param/similarity.param') ##impact??
dummy_atom1 = [0.8, 0.7, 0.9]
dummy_atom2 = [0.6, 0.5, 0.4]
dummy_atom3 = [0.1, 0.2, 0.3]
UNK_INT = 6
DUMMY_COORD = 17.21606
DUMMY_ANGLE = 45.73961
DUMMY_DIHEDRAL = 13.21566
OPLS_CONVERSION_FILE = 'param.dat'
OPLS_VERSION = '14'




def find_tors_in_log(filename):
    """
    Return the atom with torsions from a log file

    Input:
      filename: logfile

    Output:
      torsions: list of atoms with torsions 
      example--> [[1,2],[4,6],[9.10]]
      (Atom 1 and 2, 4 and 6, 9 and 10 are respectevly from the same functional group)
    """
    f = open(filename, 'r')
    out_tors = []
    read_next = 0
    while f:
        line = f.readline()
        if line == "":
            break
            # With unnatural residues, sometime we expect to have a failure
        #    a = re.search(r'[E][Rr][Rr][Oo][Rr]',line)
        #    if(a):
        #      raise Exception ("Error running Macromodel to find bonds\n"+line);
        a = re.search(r'Found Tors for atoms\s+(\d+)\s+(\d+)', line)
        if (a):
            b = [int(a.group(1)) - 1, int(a.group(2)) - 1];
            b.sort()
            out_tors.append(b)
        if (read_next == 1):
            a = re.search(r'(\d+)\s+(\d+)', line)
            read_next = 0
            if (a):
                b.append(int(a.group(1)) - 1)
                out_tors.append(b)
        a = re.search(r'adding a ring closure using atoms:\s+(\d+)\s+(\d+)', line)
        if (a):
            b = [int(a.group(2)) - 1];
            read_next = 1;
    f.close()
    return out_tors


####################################
def find_tors_in_rings(tors, maefile):
    # Find all the flexible rings in the mae file that are involved in a torsion.  As the torsions in 
    # a closed ring are heavily linked, it is necessary to include all torsions of a ring even if they
    # are double bonds.  In constrained rings we easily see changes of +/- 30 deg which has a BIG effect
    # downstream and will often prevent ring closure. Ideally we could determine which torsions change in
    # the output file, but that is for another day.   
    st1 = structure.StructureReader(maefile).next()
    cur_ring_num = 0
    out_tors = []
    out_ring = []
    for ring in st1.ring:
        ring_atoms = ring.getAtomList()
        include_ring = 0
        for t in tors:
            for atom1 in ring_atoms:
                if (t[0] + 1 == atom1):
                    for atom2 in ring_atoms:
                        if (t[1] + 1 == atom2):
                            include_ring = 1

        if (include_ring == 1):
            cur_ring_num = cur_ring_num + 1
            for atom1 in ring_atoms:
                for atom2 in ring_atoms:
                    if (atom2 > atom1 and st1.getBond(st1.atom[atom1], st1.atom[atom2]) != None):
                        #            print "ADD TORS ",st1.atom[atom1].property['s_m_pdb_atom_name'],st1.atom[atom2].property['s_m_pdb_atom_name']
                        out_tors.append([atom1 - 1, atom2 - 1])
                        out_ring.append(cur_ring_num)

    # Merge together rings which share atoms
    for i in range(len(out_tors)):
        for j in range(i + 1, len(out_tors)):
            if (out_ring[i] != out_ring[j] and
                    (out_tors[i][0] == out_tors[j][0] or
                             out_tors[i][0] == out_tors[j][1] or
                             out_tors[i][1] == out_tors[j][0] or
                             out_tors[i][1] == out_tors[j][1])):
                old_ring = max([out_ring[i], out_ring[j]])
                new_ring = min([out_ring[i], out_ring[j]])
                for k in range(len(out_tors)):
                    if out_ring[k] == old_ring: out_ring[k] = new_ring
    # Renumber rings (fill in gaps)
    if (len(out_ring) > 0):
        for ring_number in range(1, max(out_ring) + 1):
            blanks = 0
            while (not (ring_number + blanks in out_ring) and
                               ring_number + blanks < max(out_ring)): blanks = blanks + 1
            if (blanks > 0):
                for i in range(len(out_tors)):
                    if (out_ring[i] > ring_number):
                        out_ring[i] = out_ring[i] - blanks
    return out_tors, out_ring


####################################
def remove_tors(tors1, tors2):
    """
    return the tors2 from tors 1
    """
    out_tors = []
    for torsion1 in tors1:
        found = 0
        for torsion2 in tors2:
            if torsion1 == torsion2:
                found = 1
        if found == 0:
            out_tors.append(torsion1)
    return out_tors


####################################
def add_tors(tors1, tors2):
    """
    adds tors2 to tors 1
    """
    out_tors = tors1
    for torsion_i in tors2:
        found = 0
        for torsion_j in out_tors:
            if torsion_i == torsion_j:
                found = 1
        if found == 0:
            out_tors.append(torsion_i)
    return out_tors


####################################
def intersect_tors(tors1, tors2):
    # return the intersection of tors1 and tors2
    out_tors = []
    for torsion_i in tors1:
        found = 0
        for torsion_j in tors2:
            if torsion_i == torsion_j:
                found = 1
        if found == 1:
            out_tors.append(torsion_i)
    return out_tors


####################################
def mass_of_element(element):
    if (element == 'H'):
        return 1.01
    elif (element == 'D'):
        return 2.0
    elif (element == 'C'):
        return 12.01
    elif (element == 'N'):
        return 14.01
    elif (element == 'O'):
        return 16.00
    elif (element == 'F'):
        return 19.00
    elif (element == 'F'):
        return 19.00
    elif (element == 'P'):
        return 30.97
    else:
        return 10.0


####################################
def find_names_in_mae(filename, undersc=False):

    """
    :param: the file name to parse.
    :rtype: list with atom names
    This function parses a .mae file to extract the atom names from the bond section.
    This section follows the pattern:
      bond[(0-9)*]{
      keywords
      :::
      values
      }
      The first line specifies the section followed by the number of atoms enclosed with [].
      The following lines specify which value is stored in each column in the values section,
      there's one keyword (column name) by line.
      The values section has as many lines as atoms in the molecule, the number specified next
      to "bonds". Each line should have as many fields (separated by blank spaces) as keywords
      were present, and these fields should be in the same order as the keywords.
    the function ignores all the lines but the ones in the bond section.
    It reads the keywords in the bond section into a list and then looks in the list for the
    fields containing the atom_pdb_name and the pdb_residue_name, using the regular expressions:
    '.*pdb_*atom_*name' and '.*pdb_*res[idue_]*name.*', respectively, it's case insensitive, and
    they should be easy to modify
    J.M.I.F
    """
    f = open(filename, "r")
    names = []
    keywords = []
    while f:  # Find Bond Section
        line = f.readline()
        if line == "" or re.search(r'm_atom\[\d+\]', line):
            break
    while f:
        # This
        line = f.readline()
        if line == "" or re.search(':::', line):
            break
        else:
            if "index" in line:
                keyword = "index"
            else:
                keyword = line.strip()
            keywords.append(keyword)
    while f:
        # Read in Atomnames using a list with the keywords to use the right index.
        line = f.readline()
        if line == "" or re.search(':::', line):
            break
        mae_atom_values = parse_mae_line(line)
        residue_names_index = None
        atom_names_index = None
        for index, key in enumerate(keywords):
            # This block looks for the right indexes in the keywords using regular expressions,
            # if it stops working modify the expressions to match the newer and older formats. J.M.I.F
            if re.search(r'.*pdb_*res[idue_]*name.*', key, re.IGNORECASE):
                residue_names_index = index
            elif re.search (r'.*pdb_*atom_*name', key, re.IGNORECASE):
                atom_names_index = index
        if residue_names_index is None or atom_names_index is None:
            error_message = "The keywords in the atom section form the .mae file don't match the regular " \
                            "expressions currently implemented."
            raise Exception (error_message)
        if len(mae_atom_values) >= 13:
            ace = re.search('ACE', mae_atom_values[residue_names_index])  #added by mcclendon:a ligand or modified
            # residue can't have reserved name for protein capping group ACE or NMA
            nma = re.search('NMA', mae_atom_values[residue_names_index])  #added by mcclendon
            if ((not ace) and (not nma)):
                atomname = mae_atom_values[atom_names_index]
                names.append(atomname.strip())
    f.close()

    if(undersc):
      names = ['_{}_'.format(name.strip()) for name in names]

    return names
##################################################

def check_repite_names(atomnames):
  """
  Check if the mae_file contains any
  repited name. If it's like this
  raise and error.
  """
  atoms_repited = []
  for i, atomname_target in enumerate(atomnames):
    index = i
    while(index!=0):
      index-=1
      if(atomname_target == atomnames[index]):
        identifier = ' Atom:{}, AtomType:{}'.format(i, atomname_target)
        raise Exception(ERROR_ATOMTYPES + identifier)  
  

#################################################333

def replace_vdwr_from_library(rotamer_library):
  """
    Check if the rotamer library 
    contains any vdw radius= 0
    and replace them by 0.5000
  """
  found = False
  lines = []
  radius_vdw_info, start_index, end_index = parse_nonbonded(rotamer_library)
  for i, rdw_line in enumerate(radius_vdw_info):
    NBOND_info = rdw_line.split()
    rdw = float(NBOND_info[1])/2.0
    
    if(rdw == 0):
      warnings.warn("Van der Waals of atom {} = 0 changed to default 0.5".format(NBOND_info[0]))
      NBOND_info[1] = DEFAULT_RADIUS_VDW
      found = True
    lines.append(NBOND_info)

  if found:
    rvdw_change(rotamer_library, lines, start_index, end_index)



def parse_nonbonded(rotamer_library):
  """
    Find Non bonded parameters inside
    the rotamer's library
  """
  NBN_lines = []
  with open(rotamer_library, 'r') as f:
    lines = f.readlines() 
    for i, line in enumerate(lines):
      line = line.strip('\n')
      if(line == 'NBON'):
        start_index = i+1
      elif(line == 'BOND'):
        end_index = i
    try:
      for i in range(start_index, end_index):
        NBN_lines.append(lines[i])
      return NBN_lines, start_index, end_index
    except NameError:
      print(ERROR_ROTAMER_LIB)
    
def rvdw_change(rotamer_library, radius_vdw_info, start_index, end_index):
  """
    Change all radius vanderwals 0 to 0.5
    from the rotamer's library file
  """
  with open(rotamer_library, 'r') as f:
    lines = f.readlines()
    for i, new_line in enumerate(radius_vdw_info):
      lines[start_index + i] = '{0:>5} {1:>8} {2:>8} {3:>10} {4:>8} {5:>8} {6:>13} {7:>13}\n'.format(*new_line)

  with open(rotamer_library, 'w') as f:
    f.write(''.join(lines))




####################################
def parse_mae_line(line):
    ''' 
  Notice that this function is the same than MaeFileBuilder.__tokenizeLine . Try to delete this one once it will be no longer needed
  '''

    output = []
    while (len(line) > 0):
        a = re.search(r'^\s*(\S+)(.*)', line)
        if (a):
            b = re.search(r'\"', a.group(1))
            if (b):
                a = re.search(r'^\s*\"([^\"]*)\"(.*)', line)
                if (a):
                    output.append(a.group(1))
                    line = a.group(2)
                else:
                    print("Error in mae formating\n")
                    return -1;
            else:
                output.append(a.group(1))
                line = a.group(2)
        a = re.search(r'^\s*$', line)
        if (a): 
          break;
    return output


####################################
def find_mass_names(names):
    mass = []
    for i in names:
        mass.append(mass_of_element(i[1]))
    return mass


####################################
def find_bonds_in_mae(filename):
    """
    Search for bonds in MAE

    Input:
      filename: MAEfile

    Output:
      Bonds: Bonds in MAEfile
      ex--> [[0, 1], [0, 5], [0, 6]]
            (Bond between atom 0 and 1, etc...)
    """
    f = open(filename, "r")
    out_bond = []
    keywords = []
    while f:  # Find Bond Section
        line = f.readline()
        if line == "" or re.search('m_bond', line):
            break
    while f:  # Advance to numbers
        line = f.readline()
        if line == "" or (re.search(':::', line)):
            break
        else:
            if "index" in line:
                keyword = "index"
            else:
                keyword = line.strip()
            keywords.append(keyword)
    # print keywords
    while f:  # Read in Bonds
        line = f.readline()
        if line == "" or re.search(':::', line):
            break
        a = re.match(r'\s+\d+\s+(\d+)\s+(\d+)\s+(\d+)', line)
        if (a and int(a.group(3)) > 0 ):
            b = [int(a.group(1)) - 1, int(a.group(2)) - 1]
            b.sort()
            out_bond.append(b)
    f.close()
    return out_bond


####################################
def find_connected(atom, bonds, assign):
    """
    Find and assign the same "group (number)" to all the atoms connected to atom

    Input:
      atom: atom to look connections from
      bonds: list of all bonds
      assign: list of atoms with numbers assigned in order to cluster them
      ex: -->[1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 2, 2, 2, 2, 3, 3, 3]
      Atom 1 connected to 2 and 3 so they are in the same group specified by the number 1, etc...

    """
    for i in range(len(bonds)):
        if (bonds[i][0] == atom and assign[bonds[i][1]] == 0 ):
            assign[bonds[i][1]] = assign[atom]
            find_connected(bonds[i][1], bonds, assign)
        if (bonds[i][1] == atom and assign[bonds[i][0]] == 0 ):
            assign[bonds[i][0]] = assign[atom]
            find_connected(bonds[i][0], bonds, assign)


####################################
def assign_ligand_groups(tors, all_bonds, n_atoms):
    """
    Cluster atoms in groups depending whether or not they are connected

    Input:
      tors: all torsions of the ligand
      all_bonds: all ligand bonds
      n_atoms: number of atoms of the ligand

    Output:
      assign= List of Cluster of atoms depending whether or not they are connected
      ex: -->[1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 2, 2, 2, 2, 3, 3, 3]
      Atom 1 connected to 2 and 3 so they are in the same group specified by the number 1, etc...
    """

    bonds = remove_tors(all_bonds, tors)  # fixed bonds
    n_assign = 0
    c_group = 0
    assign = range(n_atoms)
    for i in range(n_atoms):
        assign[i] = 0
    done = 0
    while (done == 0):
        done = 1
        # Find the first unassigned atom and assign to a new group
        for i in range(n_atoms):
            if ( assign[i] == 0):
                c_group = c_group + 1;
                done = 0
                assign[i] = c_group
                # Find all the atoms that connect to it
                find_connected(i, bonds, assign)
                break
    return assign


####################################
def find_largest_ligand_group(assign, mass):
    max_mass = -1
    max_group = -1
    for i in range(max(assign) + 1):
        mass_group = 0
        for j in range(len(assign)):
            if (assign[j] == i):
                mass_group = mass_group + mass[j]
        if (mass_group > max_mass):
            max_mass = mass_group
            max_group = i

    output = []
    for j in range(len(assign)):
        if (assign[j] == max_group):
            output.append(j)
    return output


####################################
def convert_num_to_name(number, names):
    output = []
    for i in number:
        output.append(names[i])
    return output


####################################
def convert_name_to_num(mynames, names):
    output = []
    for thisname in mynames:
        output.append(names.index(re.sub('_', ' ', thisname)))
    return output


####################################
def assign_rank_group(atom_num, assign, rank, rank_num):
    """
    Input:
      atom_num: Atom respectively which we are producing the rank from
      assign: Group of molecules grouped on clusters dependiong on the ligand connectivity
      rank: list of numbers for each atom which will show 
            which atoms are closer to the atom we are making the rank from.
      rank_number: 0
    """

    # Add assignments for all atom in the same group as atom_num
    for i in range(len(assign)):
        if (assign[i] == assign[atom_num]):
            rank[i] = rank_num
    return rank


####################################
def min_value(array):
    if (len(array) <= 0):
        min_value = [];
        return min_value;
    min_value = array[0]
    for i in array:
        if (i < min_value):
            min_value = i
    return min_value


####################################
def max_value(array):
    if (len(array) <= 0):
        max_value = [];
        return max_value;
    max_value = array[0]
    for i in array:
        if (i > max_value):
            max_value = i
    return max_value


####################################
def assign_rank(bonds, assign, atom_num):
    """
    Define a list of ranks for each grup of atoms or cluster in assign
    which will show which atoms are closer to the group.
    As small is the number of the rank as close to the atom we are making the rank from it will be.

    Ex:
    Rank From number 1:
      Start: [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
      Middle: [0, 0, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1] (look for the closest atoms)
      Final: [0, 0, 0, 1, 2, 3, 4, 3, 2, 3, 4, 5, 5, 3, 4, 5, 4, 6, 6, 6] As small is the number as close to atom number 2 is positioned on the ligand.
    
    Input:
      bonds: Ligand connectivity
      assign: List of atoms with numbers assigned in order to cluster them
      atom_num: number of the atom to calculate the rank respect to.
    Output:
      Rank: list of numbers for each grup of atoms or cluster in assign
            which will show which atoms are closer to the group.
    """

    rank = []
    num_assign = 1;
    for i in range(len(assign)):
        rank.append(-1)
    rank = assign_rank_group(atom_num, assign, rank, 0)  # assigns atom_num a rank of zero
    while (min_value(rank) < 0):  # do until all are assigned
        cur_rank = max_value(rank)  #the loop is over ranks
        changed = 1
        while (changed == 1):  #while we haven't assigned any atoms
            changed = 0
            for i in range(len(bonds)):  #bonds is a list of pairs of atoms that are bonded
                if ( rank[bonds[i][0]] == cur_rank and rank[bonds[i][1]] < 0):
                    rank[bonds[i][1]] = rank[bonds[i][0]] + 1
                    changed = 1
                if ( rank[bonds[i][1]] == cur_rank and rank[bonds[i][0]] < 0):
                    rank[bonds[i][0]] = rank[bonds[i][1]] + 1
                    changed = 1
    return rank


####################################
def assign_group(bonds, rank):
    """
      With the core atom and its rank defined
      grouped the atoms in cluster or rotatable chains.

      Input:
        bonds: all bonds
        rank: core atom rank

      Output:
        group: groups of rotatable chains
    """
    group = []
    cur_group = -1
    for i in range(len(rank)):
        if (rank[i] == 0):
            group.append(-1)  # core atoms
        else:
            group.append(-2)  # unknown group
    while ( min_value(group) < -1 ):
        # Find an atom of rank 1 that is not assigned 
        cur_atom = -1
        for i in range(len(rank)):
            if (rank[i] == 1 and group[i] == -2):
                cur_atom = i;
                cur_group = cur_group + 1;
                group[i] = cur_group;
                break;
        atoms_in_group = 1
        if (cur_atom == -1):
          print("ERROR 2042!!")
        edit_any = 1;
        # Find all connected atoms of rank >= 1
        while (edit_any == 1):
            edit_any = 0
            for i in range(len(bonds)):
                if ( group[bonds[i][0]] == cur_group and group[bonds[i][1]] == -2):
                    group[bonds[i][1]] = cur_group;
                    atoms_in_group = atoms_in_group + 1;
                    edit_any = 1;
                if ( group[bonds[i][1]] == cur_group and group[bonds[i][0]] == -2):
                    group[bonds[i][0]] = cur_group;
                    atoms_in_group = atoms_in_group + 1;
                    edit_any = 1;
        if (atoms_in_group == 1):  # if there is one atom then it goes with the core
            for i in range(len(group)):
                if (group[i] == cur_group): group[i] = -1;
            cur_groups = cur_group - 1
    return group


####################################
def order_atoms(bonds, tors, back_tors, assign, rank, group):
    """
      Order the core atomss by connectivity. (distance)
      Update rank & group information with them.

      Example:
      order: [1,3,5,7,8,9] (starting from the core atom we ordem they by connectivity.)
      parent: [-1,1,3,1,3,5,7] (We keep track of the parents so we know the atom 3 
                                is connected to the 1, the 5 to 1, the 7 to the 3rd.)

    Input: 
      bonds: Connectivity
      tors: Atom with torsions
      back_tors: Backbone_torsions
      assign: list of numbers representing proximity.
      rank: list of numbers representing atom clustering depending on connecitvity

    Output:
      ordering: Order of atoms starting from the core by connectivity.
      out_parents: atom parent of the one at the same position inside the ordering list
      out_rank: rank of the ordering atoms
      out_group: Group of the ordering atoms


    """
    ordering = []
    parent = []
    # once an atom is added its assign value is set to -1
    # build ordering and parents of backbone (not very inteligent)
    [start_atom, num_core] = get_start_atom(tors + back_tors, rank)
    if (start_atom < 0):
        raise Exception("Core must be at least two atoms\n")
    ordering.append(start_atom)
    parent.append(-1)
    assign[start_atom] = -1
    while (len(ordering) < num_core):
        for i in range(len(bonds)):
            if (rank[bonds[i][0]] == 0 and rank[bonds[i][1]] == 0 and assign[bonds[i][0]] == -1 and assign[
                bonds[i][1]] > -1 ):
                ordering.append(bonds[i][1])
                parent.append(bonds[i][0])
                assign[bonds[i][1]] = -1
                break
            if (rank[bonds[i][1]] == 0 and rank[bonds[i][0]] == 0 and assign[bonds[i][1]] == -1 and assign[
                bonds[i][0]] > -1 ):
                ordering.append(bonds[i][0])
                parent.append(bonds[i][1])
                assign[bonds[i][0]] = -1
                break
    # Loop through the groups one at a time
    # Assign sidechain parent must be of rank exactly one above
    for grp in range(-1, max_value(group) + 1):
        cur_rank = 0
        while (len(ordering) < len(assign)):
            edit_any = 0
            for i in range(len(bonds)):
                if (rank[bonds[i][1]] - rank[bonds[i][0]] == 1 and assign[bonds[i][0]] == -1 and assign[
                    bonds[i][1]] > -1 and rank[bonds[i][0]] == cur_rank and group[bonds[i][1]] == grp):
                    ordering.append(bonds[i][1])
                    parent.append(bonds[i][0])
                    assign[bonds[i][1]] = -1
                    edit_any = 1;
                    break
                if (rank[bonds[i][0]] - rank[bonds[i][1]] == 1 and assign[bonds[i][1]] == -1 and assign[
                    bonds[i][0]] > -1 and rank[bonds[i][1]] == cur_rank and group[bonds[i][0]] == grp):
                    ordering.append(bonds[i][0])
                    parent.append(bonds[i][1])
                    assign[bonds[i][0]] = -1
                    edit_any = 1;
                    break
            if (edit_any == 0): cur_rank = cur_rank + 1
            if (cur_rank > max_value(rank)): break

    # Adjust parent list so it matches 
    out_parent = []
    out_rank = []
    out_group = []
    for i in range(len(parent)):
        out_parent.append(-100)
        out_rank.append(-100)
        out_group.append(-100)
    for i in range(len(parent)):
        out_rank[i] = rank[ordering[i]]
        out_group[i] = group[ordering[i]]
        if (parent[i] < 0):
            out_parent[i] = parent[i]
        else:
            out_parent[i] = ordering.index(parent[i])

    return ordering, out_parent, out_rank, out_group


####################################
def get_start_atom(tors, rank):
    """
    Get the number of atoms inside the core
    and the start atom of the core,
    which is the first atom of the cluster with no torsions.

    Input:
      tors: all atom with torsions
      rank: Rank: list of numbers for each grup of atoms or cluster in assign
            which will show which atoms are closer to the group.

    Output:
      start_atom: First atom of the cluster (group of atoms bonded together) 
                  that doesn not have torsions.
      num_core: number of atoms inside the cluster.

    Ex--> tors: [[1, 3], [8, 9]]
          rank: [0, 0, 0, 2, 3]
          The core has 3 atoms and the strating atom
          is the number 2 because it does not have torsions.
    """
    num_core = 0
    start_atom = -100
    for i in range(len(rank)):
        if (rank[i] == 0):
            num_core = num_core + 1
    for i in range(len(rank)):
        if (rank[i] == 0):
            in_tors = 0
            for t in tors:
                if (t[0] == i or t[1] == i):
                    in_tors = 1
            if (in_tors == 0):
                start_atom = i
                break
    return start_atom, num_core


####################################
def EliminateBackboneTors(in_tors, in_tors_ring_num, in_zmat_atoms, rank):
    tors = [];
    tors_ring_num = [];
    zmat_atoms = [];
    for i in range(len(in_tors)):
        if (rank[in_zmat_atoms[i]] > 1):
            tors.append(in_tors[i]);
            tors_ring_num.append(in_tors_ring_num[i])
            zmat_atoms.append(in_zmat_atoms[i])
    return tors, tors_ring_num, zmat_atoms


####################################
def FindCore_GetCoreAtom(tors, bonds, natoms, user_core_atom, back_tors, use_mult_lib,debug=False):
    """
    Search for the core atom wich maximes the torsions
    and it's in the center of the ligand??

    Input:
      tors: Atoms with torsions
      bonds: connectivity
      natoms: number of atoms of the ligand
      user_core_atom: predefined user core atom
      back_tors: user backbone atoms with torsion
      use_mult_lib: User decision to use or not multiples libraries.

    Ouput:
      core_atom: atom which will be the center of the core
      assign: List of atoms with numbers assigned in order to cluster them
      rank:  list of numbers for each atom which will show 
            which atoms are closer to the atom we are making the rank from.
      group: List that classifies the atoms in groups depending on connectivity
              to later on make different multiple libraries.
    """

    # Split into sections not seperated by rotatable bonds
    assign = assign_ligand_groups(tors, bonds, natoms)
    if debug:
        print(' -- ligand groups assigned.')
    if (user_core_atom > 0):
        print(' -- r')
        core_atom = user_core_atom - 1
    else:
        # print ' -- here3'
        min_rank = 10000
        core_atom = -1
        for atom_num in range(natoms):
            # print bonds, assign, atom_num
            rank = assign_rank(bonds, assign, atom_num)
            if debug:
                print(' -- rank assigned')
            [start_atom, num_core] = get_start_atom(tors + back_tors, rank)
            if (max_value(rank) < min_rank and start_atom >= 0):
                core_atom = atom_num
                min_rank = max_value(rank)
                # mass = find_mass_names(atom_names)
            #    group = find_largest_ligand_group(assign,mass) 
            #    out_names = convert_num_to_name(group,atom_names)
    if debug:
        print(' -- else finished')
    rank = assign_rank(bonds, assign, core_atom)
    if debug:
        print(' -- ranks assigned')
    if (use_mult_lib):
        group = assign_group(bonds, rank)
    else:
        group = []
        for i in range(natoms):
            group.append(0)
    return core_atom, assign, rank, group


####################################
def FindCore_GetFurthestAtom(tors, bonds, natoms, user_core_atom, back_tors, use_mult_lib):
    # Split into sections not seperated by rotatable bonds
    assign = assign_ligand_groups(tors, bonds, natoms)
    if (user_core_atom > 0):
        core_atom = user_core_atom - 1
    else:
        min_rank = 10000
        core_atom = -1
        for atom_num in range(natoms):
            rank = assign_rank(bonds, assign, atom_num)
            [start_atom, num_core] = get_start_atom(tors + back_tors, rank)
            if (max_value(rank) < min_rank and start_atom >= 0):
                core_atom = atom_num
                min_rank = max_value(rank)
                # mass = find_mass_names(atom_names)
            #    group = find_largest_ligand_group(assign,mass) 
            #    out_names = convert_num_to_name(group,atom_names)
    rank = assign_rank(bonds, assign, core_atom)
    if (use_mult_lib):
        group = assign_group(bonds, rank)
    else:
        group = []
        for i in range(natoms):
            group.append(0)
    return core_atom, assign, rank, group


####################################
def assign_bonds_to_groups(tors, group):
    """
    Make a group for each torsion bond
    and keep track of how many members

    Finally it returns the biggest group.

    Input:
      Tors: atoms with torsions
      Group: Atoms grouped by proximity

    Output:
      output: lit of group_numbers
      big_grup: biggest group
      nbig_group: members on the biggest group
    """
    output = []
    big_group = -1
    nbig_group = 0
    ngroup = max(group)
    ngroup_members = []
    ##Hauria danar un mes??
    for i in range(ngroup + 1):
        ngroup_members.append(0)
    for t in tors:
        group_number = max(group[t[0]], group[t[1]])
        output.append(group_number)
        if (group_number >= 0):
            ngroup_members[group_number] = ngroup_members[group_number] + 1
    for i in range(ngroup + 1):
        if (ngroup_members[i] > nbig_group):
            nbig_group = ngroup_members[i]
            big_group = i
    return output, big_group, nbig_group


##########################################
##########################################
# The following are for unnatural amino acids
# Chris McClendon, Jacobson Group

###################################
def FindAmideNitrogen(mae_file):
    atom_names = find_names_in_mae(mae_file)
    if (len(atom_names) <= 0):
      print("NO ATOMS FOUND IN MAE FILE")
    core_atom = -1
    for atom_num in range(len(atom_names)):
        if (atom_names[atom_num] == ' N  '):
            core_atom = atom_num
            # print "_N__ atom:"+str(core_atom)
    if (core_atom == -1): raise Exception("cannot find Amide Nitrogen _N__ in maefile\n")

    return core_atom


###################################
def FindAmideHydrogen(mae_file):
    atom_names = find_names_in_mae(mae_file)
    if (len(atom_names) <= 0): 
      print("NO ATOMS FOUND IN MAE FILE")
    nh_atom = -1
    for atom_num in range(len(atom_names)):
        if (atom_names[atom_num] == ' H  '): nh_atom = atom_num
    if (nh_atom == -1): raise Exception("cannot find Amide Nitrogen _N__ in maefile\n")
    # print "_H__ atom:"+str(nh_atom)
    return nh_atom


###################################
def FindHCAlpha(mae_file, raise_exception=True):
    atom_names = find_names_in_mae(mae_file)
    if (len(atom_names) <= 0):
      print("NO ATOMS FOUND IN MAE FILE")
    ha_atom = -1
    for atom_num in range(len(atom_names)):
        if (atom_names[atom_num] == ' HA '): ha_atom = atom_num
    if (ha_atom == -1):
        if (raise_exception == True): raise Exception("cannot find HCalpha  _HA_ in maefile\n")
    # print "_HA_ atom:"+str(ha_atom)
    return ha_atom


###################################
def FindCAlpha(mae_file, raise_exception=True):
    atom_names = find_names_in_mae(mae_file)
    if (len(atom_names) <= 0): 
      print("NO ATOMS FOUND IN MAE FILE")
    bb_atom = -1
    for atom_num in range(len(atom_names)):
        if (atom_names[atom_num] == ' CA '): bb_atom = atom_num
    if (bb_atom == -1):
        if (raise_exception == True): raise Exception("cannot find Calpha  _CA_ in maefile\n")
    # print "_CA_ atom:"+str(bb_atom)
    return bb_atom


###################################
def FindC(mae_file, raise_exception=True):
    atom_names = find_names_in_mae(mae_file)
    if (len(atom_names) <= 0):
      print("NO ATOMS FOUND IN MAE FILE")
    bb_atom = -1
    for atom_num in range(len(atom_names)):
        if (atom_names[atom_num] == ' C  '): bb_atom = atom_num
    if (bb_atom == -1):
        if (raise_exception == True): 
          raise Exception("cannot find _C__  in maefile\n")
    # print "_C__ atom:"+str(bb_atom)
    return bb_atom


###################################
def FindO(mae_file):
    atom_names = find_names_in_mae(mae_file)
    if (len(atom_names) <= 0):
      print("NO ATOMS FOUND IN MAE FILE")
    bb_atom = -1
    for atom_num in range(len(atom_names)):
        if (atom_names[atom_num] == ' O  '): bb_atom = atom_num
    if (bb_atom == -1): raise Exception("cannot find _O__  in maefile\n")
    # print "_O__ atom:"+str(bb_atom)
    return bb_atom


###################################
def Buildup_Connected(bonds, length, connected):  # takes in dictionary object connected
    for i in range(length):
        connected[i] = []
        for j in range(len(bonds)):
            if (bonds[j][0] == i and (not (bonds[j][1] in connected[i]))): connected[i].append(bonds[j][1])
            if (bonds[j][1] == i and (not (bonds[j][0] in connected[i]))): connected[i].append(bonds[j][0])
    return connected


###################################
def Order_Atoms_AA(bonds, tors, assign, rank, group, mae_file):
    # Chris McClendon, Jacobson Group
    ordering = []
    parent = []
    # assumes at least one atom will be a core atom with rank == 0
    # once an atom is added its assign value is set to -1
    # build ordering and parents of backbone (not very inteligent)
    num_core = 0
    start_atom = -100
    for i in range(len(rank)):
        if (rank[i] == 0):
            num_core = num_core + 1
    start_atom = -1
    Calpha_atom = -1
    HCAlpha_atom = -1
    C_atom = -1
    start_atom = FindAmideNitrogen(mae_file)
    Calpha_atom = FindCAlpha(mae_file)
    HCalpha_atom = FindHCAlpha(mae_file)
    C_atom = FindC(mae_file)
    backbone = (start_atom, Calpha_atom, C_atom)
    for stuff in backbone:
        if (stuff == -1): raise Exception("cannot identify N, Calpha, or C backbone atoms!")

    def OnlyParentConnected(atom, current_parent, connected):
        myval = 0
        if (connected.has_key(atom) and connected[atom][0] == current_parent):
            if (len(connected[atom]) == 1):
                myval = 1
        return myval

    def TryAddAtom(atom, current_parent, assign):
        #print("trying to add atom:"+str(atom))
        if (assign[atom] != -1):  #if atom not already assigned
            ordering.append(atom)
            parent.append(current_parent)
            assign[atom] = -1


        #have to modify to handle Pro analogs


    def AppendSideChains(current_parent):
        #print("appending to atom:"+str(current_parent))
        #print("assign states:"+str(assign))
        assert (assign[
                    current_parent] == -1)  #loop precondition: parent is already assigned and is the current value of current_parent
        if (not (connected.has_key(current_parent))): return  #failsafe base case
        for atom in connected[current_parent]:
            if (assign[atom] > -1):
                if (OnlyParentConnected(atom, current_parent,
                                        connected)):  #handle hydrogens and other single-atom groups as the base case
                    TryAddAtom(atom, current_parent, assign)
        for atom in connected[current_parent]:  #now look for recursion
            if (assign[atom] > -1):  #if already assigned, don't go down that path
                if (atom not in backbone):  #now handle recursive case
                    TryAddAtom(atom, current_parent, assign)  #append the atom noting the current parent
                    AppendSideChains(atom)  #and then use this atom as the parent recursively

    connected = {}


    #Construct array of lists of connections
    connected = Buildup_Connected(bonds, len(rank), connected)

    #Start off with Nitrogen

    TryAddAtom(start_atom, -1, assign)
    #print("ordering before append to nitrogen\n"+str(ordering))
    #print("parent before append to nitrogen\n"+str(parent))
    AppendSideChains(start_atom)  #recurisve to include HN(s) or other atoms attached to N

    TryAddAtom(Calpha_atom, start_atom, assign)
    #print("ordering before append to Calpha\n"+str(ordering))
    #print("parent before append to Calpha\n"+str(parent))
    AppendSideChains(Calpha_atom)  #recursive to include sidechain

    TryAddAtom(C_atom, Calpha_atom, assign)
    #print("ordering before append to C \n"+str(ordering))
    #print("parent before append to C\n"+str(parent))
    AppendSideChains(C_atom)  #recursive to include O or ther modifications like OXT
    #print("ordering after append to C \n"+str(ordering))


    #  Adjust parent list so it matches
    out_parent = []
    out_rank = []
    out_group = []
    for i in range(len(parent)):
        out_parent.append(-100)
        out_rank.append(-100)
        out_group.append(-100)
    for i in range(len(parent)):
        out_rank[i] = rank[ordering[i]]
        out_group[i] = group[ordering[i]]
        if (parent[i] < 0):
            out_parent[i] = parent[i]
        else:
            out_parent[i] = ordering.index(parent[i])

    return ordering, out_parent, out_rank, out_group


####################################
def Reorder_Amide_Nitrogen_Hydrogen(ordering, parent, rank, group, nitrogen_atom, hydrogen_atom):
    # Chris McClendon, Jacobson Group
    out_parent = []
    out_rank = []
    out_group = []
    new_ordering = []
    for i in range(len(ordering)):
        new_ordering.append(-1)
    for i in range(len(new_ordering)):
        if (i > 1): new_ordering[i] = ordering[i]
        if (i == 0):
            if (ordering[0] == hydrogen_atom and ordering[1] == nitrogen_atom):
                print("Swapping Nitrogen and Hydrogen in Ordering")
                new_ordering[0] = nitrogen_atom
                new_ordering[1] = hydrogen_atom
                #parent[nitrogen_atom]=-1; parent[hydrogen_atom]=nitrogen_atom
            else:
                new_ordering[0] = ordering[0]
                new_ordering[1] = ordering[1]

    for i in range(len(parent)):
        out_parent.append(-100)
        out_rank.append(-100)
        out_group.append(-100)
    for i in range(len(parent)):
        # if(ordering[0]==hydrogen_atom):
        #   out_rank[i] = rank[new_ordering[i]]+1
        # else:
        out_rank[i] = int(rank[ordering[i]])
        out_group[i] = group[ordering[i]]
        if (parent[i] < 0):
            out_parent[i] = parent[i]
        else:
            out_parent[i] = parent[i]

    return new_ordering, out_parent, rank, out_group


#######################################
def FindCoreAA(mae_file, user_fixed_bonds, use_rings, residue_name, use_mult_lib, user_core_atom, user_tors):
    # adapted from FindCore for unnatural amino acids by Chris McClendon, Jacobson Group
    if (user_tors == []):
        tors = get_torsions_from_mae(mae_file, residue_name)
        ########CHANGE     SCHRODINGER######
        [ring_tors, ring_num] = find_tors_in_rings(tors, mae_file)
        ########CHANGE     SCHRODINGER######
        tors = remove_tors(tors, user_fixed_bonds)

        #########CHANGE     SCHRODINGER######
        ring_tors = remove_tors(ring_tors, user_fixed_bonds)
        if (len(ring_tors) < 1 ): 
          use_rings = 0
        if (use_rings == 0):
          tors = remove_tors(tors, ring_tors)
    else:
        tors = user_tors
        ##########CHANGE     SCHRODINGER######

    bonds = find_bonds_in_mae(mae_file)
    atom_names = find_names_in_mae(mae_file)
    if (len(atom_names) <= 0): 
      print("NO ATOMS FOUND IN MAE FILE")
    print("bonds\n")
    print(bonds)
    print("tors\n")
    print(tors)
    assign = assign_ligand_groups(tors, bonds, len(atom_names))
    if (user_core_atom > 0):
        core_atom = user_core_atom - 1
    else:
        core_atom = FindAmideNitrogen(mae_file)
    start_atom = FindAmideNitrogen(mae_file)
    Calpha_atom = FindCAlpha(mae_file)
    HCalpha_atom = FindHCAlpha(mae_file)
    C_atom = FindC(mae_file)
    O_atom = FindO(mae_file)
    #assign Ca, HC, and C atoms to same group as backbone nitrogen or core atom so these get designation ' M' are treated as backbone in PLOP
    assign[Calpha_atom] = assign[core_atom]
    assign[HCalpha_atom] = assign[core_atom]
    assign[C_atom] = assign[core_atom]
    assign[O_atom] = assign[core_atom]
    #now set ranks according to assigned groups
    rank = assign_rank(bonds, assign, core_atom)
    if (use_mult_lib):
        group = assign_group(bonds, rank)
    else:
        group = []
        for i in range(len(atom_names)):
            group.append(0)
    line = 'Core Atoms: '
    for i in range(len(atom_names)):
        if (rank[i] == 0):
            line = line + atom_names[i].rjust(5)
    print('\n')
    print(line)
    #  for i in range(len(atom_names)):
    #    print atom_names[i].rjust(5),group[i],rank[i]

    #  print "Core atom",core_atom,rank
    if (use_mult_lib == 1):
        print('\n')
        print('Number of groups {}'.format(max_value(group) + 1))
        for grp in range(max_value(group) + 1):
            line = 'Group %2d' % grp
            for i in range(len(atom_names)):
                if (group[i] == grp):
                    line = line + atom_names[i].rjust(5)
            print(line)

    [old_num, parent, rank, group] = Order_Atoms_AA(bonds, tors, assign, rank, group, mae_file)
    #  for i in range(len(atom_names)):
    #    print atom_names[i].rjust(5),group[i],rank[i],parent[i]


    # Assign ring numbers to torsions (0 for none) (backbone ring torsions not assigned)
    tors_ring_num = []
    for t in tors: 
        tors_ring_num.append(0)
    if (use_rings == 1):
        for i in range(len(tors)):
            for j in range(len(ring_tors)):
                if (tors[i][0] == ring_tors[j][0] and tors[i][1] == ring_tors[j][1]):
                    tors_ring_num[i] = ring_num[j]

    return old_num, parent, rank, tors, use_rings, group, tors_ring_num


####################################
def ReorderTorsionsAA(tors, ordering):  # by Chris McClendon, Jacobson Group
    oldtors = []
    newtors = []
    inv_ordering = []
    for i in range(len(tors)): oldtors.append(tors[i])  #copies by value

    #Assuming atom ordering is PLOP tree.F compatible, we simply go down the atom list andd write down the atoms having a torsion

    #Every torsion represents a real rotatable bond
    #We need to know whether the first or second atom of a torsion is in the direction of the tree starting from the amide nitrogen
    #The atom with is further down in the ordering will represent the direction of the tree

    #create an inverse ordering map that takes a mae_atom_num and gives it's location in the ordered list
    for i in range(len(ordering)):
        inv_ordering.append(0)
    for i in range(len(ordering)):
        inv_ordering[ordering[i]] = i

    #now reorder the torsions in a PLOP tree.F compatible way to prepare for TetherRotBonds.output_rotbonds()
    for i in range(len(ordering)):
        if (i > 0):  #skip amide nitrogen, which we assume is always first, even for a peptoid
            foundtors = 0
            #print "oldtors:" + str(oldtors)
            for mytors in oldtors:
                if (inv_ordering[mytors[1]] > inv_ordering[mytors[0]]):  # then tree direction is mytors[0] -> mytors[1]
                    if (mytors[1] == ordering[i]):
                        newtors.append([mytors[0], ordering[i]])
                        oldtors = remove_tors(oldtors, [mytors])
                        foundtors == 1
                        break
                if (inv_ordering[mytors[1]] < inv_ordering[mytors[0]]):  # then tree direction is mytors[1] -> mytors[0]
                    if (mytors[0] == ordering[i]):
                        newtors.append([mytors[1], ordering[i]])
                        oldtors = remove_tors(oldtors, [mytors])
                        foundtors == 1
                        break
    print('Found the core')
    return newtors


####################################################
# End Unnatural Amino Acid Code
####################################################
def FindCore(mae_file, user_fixed_bonds, use_rings, residue_name,
             use_mult_lib, user_core_atom, user_tors, back_tors, max_tors, R_group_root_atom_name):
    """

    Find a core group of atoms and cluster the others.
    
    Input:
      mae_file: ligand topology
      user_fixed_bonds: self-explanatory
      logfile: file with the connectivity and torsions
      use_rings: option to search for ring torsions
      use_mult_lib: use mult library for each cluster
      user_core_atom: create a cluster starting on these atom
      user_tors: torsions specified by the user
      back_torsions: backbone torsions specified By the user
      max_tors: maximum number of torsions in each group
      R_group_root_atom_name: ???

    Output:
      old_num: Order of atoms by connectivuity starting from the centreal atom of the core
      parent: atom connected with the one on the same position inside the old_num list
      rank: atoms in old_num ranked by proximity
      tors: Atoms with torsions (backbone excluded)
      use_rings: option to search for ring torsions
      back_tors: Atoms of the backbone with torsions
      tors_ring_num: Atoms coming from rings with torsions. (0 for none)

    example--> 
    old_num: [4, 3, 8, 5, 13, 6, 14, 7, 15, 16, 1, 0, 2, 9, 10, 11, 12, 17, 18, 19]
    parent: [-1, 0, 1, 0, 0, 3, 3, 5, 5, 7, 1, 10, 10, 2, 13, 14, 14, 16, 16, 16]
    rank: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 2, 1, 2, 3, 3, 4, 4, 4]

    - Atom 4 is close to 3 close to 8 to 5, etc...
    - Atom 3 is connected to 0, 8 to 1, 5 to 0, 13 to 0, etc...
    - Atom 4,3,8,5,13,6,14,7,15 are from the rank 0 definides as core. All other ranks are sidechains.

    """
    if (user_tors == []):
        # tors = find_tors_in_log(logfile)
        tors = get_torsions_from_mae(mae_file, residue_name)
        print(tors)
        print(" - torsions found.")

    ########CHANGE     SCHRODINGER######
        [ring_tors, ring_num] = find_tors_in_rings(tors, mae_file)
        print(" - ring torsions found.")
    ########CHANGE     SCHRODINGER######
        tors = remove_tors(tors, user_fixed_bonds)
        tors = remove_tors(tors, back_tors)
        ring_tors = remove_tors(ring_tors, user_fixed_bonds)
        if (len(ring_tors) < 1 ): use_rings = 0
        if (use_rings == 0):
            tors = remove_tors(tors, ring_tors)
        else:
            tors = add_tors(tors, ring_tors)

    else:
        tors = user_tors
    if (max_tors != -1 and use_rings == 1):
        print("WARNING: May freeze some torsions in rings with unpredictable results")
    #  print "tors"
    #  for t in tors:
    #    print t
    #  print "ring tors"
    #  for t in ring_tors:
    #    print t
    bonds = find_bonds_in_mae(mae_file) 
    print(' -bonds found in mae file')
    atom_names = find_names_in_mae(mae_file)
    print(atom_names) 
    print(' -names found in mae file')
    if (len(atom_names) <= 0):
      print("NO ATOMS FOUND IN MAE FILE")

    if (user_core_atom == -2):  #then find the atom furthest from the R group
        R_group_root_atom_num = convert_name_to_num([R_group_root_atom_name], atom_names)
        [core_atom, assign, rank, group] = FindCore_GetCoreAtom(tors, bonds, len(atom_names), R_group_root_atom_num[0],
                                                                back_tors, use_mult_lib)
        user_core_atom = rank.index(max(rank))
 
    [core_atom, assign, rank, group] = FindCore_GetCoreAtom(tors, bonds, len(atom_names), user_core_atom, back_tors, use_mult_lib)
    print(' - core atoms found')
    [group_tors, big_group, nbig_group] = assign_bonds_to_groups(tors, group)
    print(' -bonds')
    ##############Canvir aqui torsions##########################
    while ( nbig_group > max_tors and max_tors >= 0 ):
        # Find Lowest Rank Torsion in Biggest Group
        chosen = -1
        lowest = 10000
        for i in range(len(tors)):
            if group_tors[i] == big_group:
                t = tors[i]
                this_rank = max(rank[t[0]], rank[t[1]])
                if ( this_rank < lowest ):
                    lowest = this_rank
                    chosen = i
        if (chosen < 0):
            raise Exception("ERROR Removing Torsion")
        t = tors[chosen]
        #print("Max Torsions in a sidechain group is ", nbig_group, ", Assigning Torsion to Backbone: ", atom_names[
         #   t[0]], "-", atom_names[t[1]]
        tors = remove_tors(tors, [t])
        back_tors.append(t)
        [core_atom, assign, rank, group] = FindCore_GetCoreAtom(tors, bonds, len(atom_names), user_core_atom, back_tors, use_mult_lib)
        [group_tors, big_group, nbig_group] = assign_bonds_to_groups(tors, group)


    #  print "Max Torsions is ",nbig_group  
    #  for i in range(len(tors)):
    #     t = tors[i]
    #     print atom_names[t[0]],"-",atom_names[t[1]],group_tors[i]   

    line = ' -Core Atoms (including backbone):'
    for i in range(len(atom_names)):
        if (rank[i] == 0):
            line = line + ' ' + atom_names[i]
    print(line)
    #  for i in range(len(atom_names)):
    #    print atom_names[i].rjust(5),group[i],rank[i]

    #  print "Core atom",core_atom,rank
    if (use_mult_lib == 1):
        print(' -Number of groups {}:'.format(max_value(group) + 1))
        for grp in range(max_value(group) + 1):
            line = '  -Group %2d' % (grp + 1)
            for i in range(len(atom_names)):
                if (group[i] == grp):
                    line = line + ' ' + atom_names[i]
            print(line)

    [old_num, parent, rank, group] = order_atoms(bonds, tors, back_tors, assign, rank, group)

    if (use_rings == 2): 
      tors = intersect_tors(tors, ring_tors)
      use_rings = 1

    # Assign ring numbers to torsions (0 for none) (backbone ring torsions not assigned)
    tors_ring_num = []
    for t in tors: 
      tors_ring_num.append(0)
    if (use_rings == 1):
        for i in range(len(tors)):
            for j in range(len(ring_tors)):
                if (tors[i][0] == ring_tors[j][0] and tors[i][1] == ring_tors[j][1]):
                    tors_ring_num[i] = ring_num[j]
                #  for i in range(len(atom_names)):
                #    print atom_names[i].rjust(5),group[i],rank[i],parent[i]
    return old_num, parent, rank, tors, use_rings, group, back_tors, tors_ring_num


def get_torsions_from_mae(mae_file, residue_name):
  pdb_file = residue_name + ".pdb"
  struct = structure.StructureReader(mae_file).next()
  struct.write(pdb_file)
  mol = Chem.MolFromPDBFile(pdb_file)
  torsions =  TorsionFingerprints._getBondsForTorsions(mol, True)
  torsions = [[tor[0], tor[1]] for tor in torsions]
  try:
    os.remove(pdb_file)
  except OSError:
    print("Error when calculating torsions. Be carefull not to have a {} in your current directory".format(pdb_file))
  return torsions


####################################
def ReorderTemplate(ordering, new_parent, rank, in_file, out_file, R_group_root_atom_name='None'):
    """
    Create the ain template by parsing the ain.hetgrp_ffgen
    calculating the zmatrix, ordering the connectivity,
    the bonded, non bonded & Improper rotations.

    Output: 
        names: atom names from the template file
        mae file: template file of the ligand & input for PELE.
    """

    [old_parent, zmat, temp_names] = read_zmat_template(in_file)
    #  Set up dummy values for cart add adjust parent list,ordering to match
    #  cart=[[0.8, 0.7, 0.9]]
    #  cart.append([0.6, 0.5, 0.4])
    #  cart.append([0.1, 0.2, 0.3])
    #  for i in range(len(old_parent)):
    #    old_parent[i]=old_parent[i]+3
    #    new_parent[i]=new_parent[i]+3
    #    ordering[i]=ordering[i]+3
    #  old_parent.insert(0,1);old_parent.insert(0,0);old_parent.insert(0,-1)
    #  new_parent.insert(0,1);new_parent.insert(0,0);new_parent.insert(0,-1)
    #  temp_names.insert(0,'DUM3');temp_names.insert(0,'DUM2');temp_names.insert(0,'DUM1');
    #  ordering.insert(0,2);ordering.insert(0,1);ordering.insert(0,0)
    #  zmat.insert(0,[0,0,0]);zmat.insert(0,[0,0,0]);zmat.insert(0,[0,0,0])
    #Convert to cart
    cart = int2xyz(zmat, old_parent)
    #  print 'cartesian(xyz)'
    #  for i in range(len(old_parent)):
    #     print 'X',cart[i][0],cart[i][1],cart[i][2]

    #Convert back to zmat with new parent list
    zmat = xyz2int(cart, ordering, new_parent)
    #Remove dummy atoms
    #for i in range(len(old_parent)):
    #  old_parent[i]=old_parent[i]-3
    #  new_parent[i]=new_parent[i]-3
    #  ordering[i]=ordering[i]-3
    #old_parent.pop(0);old_parent.pop(0);old_parent.pop(0)
    #new_parent.pop(0);new_parent.pop(0);new_parent.pop(0)
    #ordering.pop(0);ordering.pop(0);ordering.pop(0)
    #zmat.pop(0);zmat.pop(0);zmat.pop(0)
    #temp_names.pop(0);temp_names.pop(0);temp_names.pop(0);

    #  for i in range(len(new_parent)):
    #    print zmat[i][0],zmat[i][1],zmat[i][2]

    #Prep read in/out
    fin = open(in_file, "r")
    fout = open(out_file, "w")
    while fin:  #Get past coments
        line = fin.readline()
        fout.write(line)
        if (not re.search('^\*', line)): break
        
        

    a = re.search('^\s*\S+\s+(\d+)\s+(\d+)\s+(\d+)', line)
    if (a):
        num_atoms = int(a.group(1))
        num_bonds = int(a.group(2))
        num_theta = int(a.group(3))
    else:
        raise Exception("ERROR reading template file");
    #Read in and read out the atoms section of the template file
    at = [];
    name = [];
    mat = [];
    for i in range(num_atoms):
        line = fin.readline()
        a = re.search('^\s*(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+([\d\.\-]+)\s+([\d\.\-]+)\s+([\d\-\.]+)',
                      line)
        at.append(a.group(4));
        name.append((a.group(5)).strip('_'));
        mat.append(a.group(6));

    #Keep as mainchain the z-matrix atoms until R_group_root_atom_name,
    #the root atom of the group to be sampled, is reached
    still_mainchain = 0
    if (R_group_root_atom_name != 'None'): still_mainchain = 1

    for i in range(num_atoms):
        j = ordering[i]
        #print "Atom: "+str(i)+" Name of atom "+name[j]+" root atom name "+str(R_group_root_atom_name)+"\n"
        if (rank[i] == 0):
            rank_str = ' M'
        else:
            if still_mainchain == 1:
                if name[j] == R_group_root_atom_name:  #if we've reached the root of
                    still_mainchain = 0  #the R group, start sampling
                    rank_str = ' S'
                else:
                    rank_str = ' M'
            else:
                rank_str = ' S'
        outline = str(i + 1).rjust(5) + str(old_parent[i] + 1).rjust(6) + rank_str + '   ' + (at[j]).ljust(5) + (
        name[j]).ljust(4) + mat[j].rjust(6) + '%12.6f' % zmat[i][0] + '%12.6f' % zmat[i][1] + '%12.6f' % zmat[i][
            2] + '\n'
        fout.write(outline)
    """
    # Move to exclude table
    exclude = []
    counts = []
    while (len(counts) < len(ordering)):
        line = fin.readline();
        a = re.search('\s*(.*)', line);
        line = a.group(1)
        print(lines)
        p = re.compile(r'\W+');
        temp = p.split(line);
        counts = counts + temp
        print(counts)
    for i in range(len(ordering)):
        exclude.append([])
    for i in range(len(ordering)):
        print('ordering')
        print(len(ordering))
        line = fin.readline()
        print(line)
        a = re.search('\s*(.*)', line);
        line = a.group(1)
        print(line)
        p = re.compile(r'[\W$\s]+');
        temp = p.split(line);
        print(temp)
        iatom = conv_at(ordering, i + 1)
        print(iatom)
        for atom in temp:
            if (atom == ''): break;
            if (int(atom) == 0 ):
                pass
            else:
                #          print "Converting atom ",atom
                #          print "Ordering atom ",ordering
                jatom = conv_at(ordering, atom)
                if (int(iatom) < int(jatom) ):
                    exclude[iatom - 1].append(jatom)
                else:
                    exclude[jatom - 1].append(iatom)
    for i in range(len(ordering)):
        if (len(exclude[i]) == 0): exclude[i].append(0)
    # write out counts
    line = ''
    for i in range(len(ordering)):
        line = line + str(len(exclude[i])).rjust(4)
        if (i != 0 and (math.modf(1.0 / 16.0 * i))[0] < 0.001 and i < len(ordering)):
            line = line + '\n'
    line = line + '\n';
    fout.write(line)
    # write out exclude list
    for i in range(len(ordering)):
        line = ''
        for atom in exclude[i]:
            line = line + str(atom).rjust(6)
        line = line + '\n';
        fout.write(line)
    """
    #NONBON Region
    line = fin.readline()
    if (not re.search('NBON', line)):
        raise Exception("ERROR IN TEMPLATE FORMAT\n" + line)
    fout.write(line)
    nonbon = []
    for i in range(len(ordering)):
        line = fin.readline();
        nonbon.append(line)
    for i in range(len(ordering)):
        line = nonbon[ordering[i]]
        a = re.search('^\s*\d+(.*)', line)
        line = str(i + 1).rjust(6) + a.group(1) + '\n'
        fout.write(line)
    #BOND Region
    line = fin.readline()
    if (not re.search('BOND', line)):
        raise Exception("ERROR IN TEMPLATE FORMAT\n" + line)
    fout.write(line)
    bond = []
    for i in range(num_bonds):
        line = fin.readline();
        bond.append(line)
    for line in bond:
        a = re.search('^\s*(\d+)\s+(\d+)(.*)', line)
        line = str(ordering.index(int(a.group(1)) - 1) + 1).rjust(6) + str(
            ordering.index(int(a.group(2)) - 1) + 1).rjust(6) + a.group(3) + '\n'
        fout.write(line)
    #Theta Region
    line = fin.readline()
    if (not re.search('THET', line)):
        raise Exception("ERROR IN TEMPLATE FORMAT\n" + line)
    fout.write(line)
    theta = []
    for i in range(num_theta):
        line = fin.readline();
        theta.append(line)
    for line in theta:
        a = re.search('^\s*(\d+)\s+(\d+)\s+(\d+)(.*)', line)
        line = str(ordering.index(int(a.group(1)) - 1) + 1).rjust(6) + str(
            ordering.index(int(a.group(2)) - 1) + 1).rjust(6) + str(ordering.index(int(a.group(3)) - 1) + 1).rjust(
            6) + a.group(4) + '\n'
        fout.write(line)

    #PHI/IPHI Region
    line = fin.readline()
    if (not re.search('^PHI', line)):
        raise Exception("ERROR IN TEMPLATE FORMAT\n" + line)
    fout.write(line)
    while 1:
        line = fin.readline()
        if (not (line)): break
        if (re.search('END', line)): fout.write(line);break;
        if (re.search('IPHI', line)):
            fout.write(line)
        else:
            a = re.search('^\s*([\-\d]+)\s+([\-\d]+)\s+([\-\d]+)\s+([\-\d]+)(.*)', line)
            line = str(conv_at(ordering, a.group(1))).rjust(6) + str(conv_at(ordering, a.group(2))).rjust(6) + str(
                conv_at(ordering, a.group(3))).rjust(6) + str(conv_at(ordering, a.group(4))).rjust(6) + a.group(
                5) + '\n'
            fout.write(line)

    fin.close()
    fout.close()
    # update names
    names = []
    for i in range(len(temp_names)):
        names.append(temp_names[ordering[i]])
    return names


####################################
def conv_at(ordering, input):
    input = int(input)
    if (input < 0):
        input = -1 * (ordering.index((-1 * input) - 1) + 1)
    else:
        input = (ordering.index((   input) - 1) + 1)
    return input


####################################


def xyz2int(in_cart, in_ordering, in_parent):
    #  Set up dummy values for cart add adjust parent list,ordering to match
    cart = [dummy_atom1]
    cart.append(dummy_atom2)
    cart.append(dummy_atom3)
    parent = []
    ordering = []
    parent.insert(0, 1);
    parent.insert(0, 0);
    parent.insert(0, -1)
    ordering.insert(0, 2);
    ordering.insert(0, 1);
    ordering.insert(0, 0)
    for i in range(len(in_cart)):
        cart.append(in_cart[i])
    for i in range(len(in_parent)):
        parent.append(in_parent[i] + 3)
        ordering.append(in_ordering[i] + 3)
    #  temp_names.insert(0,'DUM3');temp_names.insert(0,'DUM2');temp_names.insert(0,'DUM1');
    zmat = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    #  for i in range(len(parent)):
    #    print cart[i][0],cart[i][1],cart[i][2]
    for i in range(3, len(parent)):
        iatom = i
        jatom = parent[iatom]
        katom = parent[jatom]
        latom = parent[katom]
        #     print "Parents(nn) ",iatom-2,jatom-2,katom-2,latom-2
        iatom = ordering[iatom]
        jatom = ordering[jatom]
        katom = ordering[katom]
        latom = ordering[latom]
        #     print "Parents(on) ",iatom-2,jatom-2,katom-2,latom-2
        xi = cart[iatom][0];
        yi = cart[iatom][1];
        zi = cart[iatom][2]
        xj = cart[jatom][0];
        yj = cart[jatom][1];
        zj = cart[jatom][2]
        xk = cart[katom][0];
        yk = cart[katom][1];
        zk = cart[katom][2]
        xl = cart[latom][0];
        yl = cart[latom][1];
        zl = cart[latom][2]
        #     print "xyz",xi,yi,zi,xj,yj,zj
        dx = xi - xj;
        dy = yi - yj;
        dz = zi - zj
        rij = math.sqrt(dx * dx + dy * dy + dz * dz)
        theta = bangle(xi, yi, zi, xj, yj, zj, xk, yk, zk)
        phi = calc_tors(xi, yi, zi, xj, yj, zj, xk, yk, zk, xl, yl, zl)
        zmat.append([rij, theta * 180.0 / math.pi, phi * 180.0 / math.pi])
    zmat = zmat[3:len(zmat)]
    return zmat


####################################
def bangle(xj, yj, zj, xi, yi, zi, xk, yk, zk):
    dijx = xj - xi
    dijy = yj - yi
    dijz = zj - zi
    dikx = xk - xi
    diky = yk - yi
    dikz = zk - zi
    vdot = dijx * dikx + dijy * diky + dijz * dikz

    if (math.fabs(vdot) > 0):
        rij = math.sqrt(dijx * dijx + dijy * dijy + dijz * dijz)
        rik = math.sqrt(dikx * dikx + diky * diky + dikz * dikz)
        xang = vdot / (rij * rik)
    else:
        raise Exception("ERROR:  zero angle in bangle");
    if (xang - 1.0 > -0.0000000001): return 0.0
    if (xang + 1.0 < 0.00000000001): return 3.13159
    return math.acos(xang)


####################################
def calc_tors(xi, yi, zi, xj, yj, zj, xk, yk, zk, xl, yl, zl):  #ported from plop "tors" in struct.F
    dijx = xi - xj
    dijy = yi - yj
    dijz = zi - zj
    dkjx = xk - xj
    dkjy = yk - yj
    dkjz = zk - zj
    dklx = xk - xl
    dkly = yk - yl
    dklz = zk - zl

    ax = dijy * dkjz - dijz * dkjy
    ay = dijz * dkjx - dijx * dkjz
    az = dijx * dkjy - dijy * dkjx
    cx = dkjy * dklz - dkjz * dkly
    cy = dkjz * dklx - dkjx * dklz
    cz = dkjx * dkly - dkjy * dklx

    rac = ax * cx + ay * cy + az * cz
    ra = ax * ax + ay * ay + az * az
    rc = cx * cx + cy * cy + cz * cz
    cosang = rac / math.sqrt(ra * rc)
    if (cosang - 1.0 > -0.00000000001):
        phi = 0.0
    elif (cosang + 1.0 < 0.00000000001):
        phi = math.pi
    else:
        phi = math.acos(cosang)

    s = dijx * cx + dijy * cy + dijz * cz
    if (s < 0): phi = - phi  # to account for phi between pi and 2pi
    return phi


####################################
def int2xyz(in_zmat, in_parent):
    #  Set up dummy values for cart add adjust parent list,ordering to match
    cart = [dummy_atom1]
    cart.append(dummy_atom2)
    cart.append(dummy_atom3)
    zmat = [[0, 0, 0]];
    zmat.append([0, 0, 0]);
    zmat.append([0, 0, 0]);
    for i in range(len(in_zmat)):
        zmat.append(in_zmat[i])
    parent = []
    for i in range(len(in_parent)):
        parent.append(in_parent[i] + 3)
    parent.insert(0, 1);
    parent.insert(0, 0);
    parent.insert(0, -1)
    # does not require that a parent atom come before the child
    # but much more efficient if it does
    epsilon = 0.000000001
    calc = []  # whether or not we have calculated a cart coordinate
    for i in range(0, 3): calc.append(1)  # we know dummies
    for i in range(3, len(zmat)): calc.append(0);cart.append([0.0, 0.0, 0.0]);
    while (min_value(calc) < 0.1):
        for i in range(3, len(zmat)):  # skip over 3 dummy atoms
            iatom = i
            jatom = parent[iatom]
            katom = parent[jatom]
            latom = parent[katom]
            if (calc[iatom] == 1): continue
            if (calc[jatom] == 0 or calc[katom] == 0 or calc[jatom] == 0): continue
            #       print "Zmat ",iatom,jatom,katom,latom
            calc[iatom] = 1  #signal we have calculated this
            rcd = zmat[iatom][0]
            thbcd = zmat[iatom][1] * math.pi / 180.0
            phabcd = 2 * math.pi - zmat[iatom][2] * math.pi / 180.0
            #       print "Zmat ",iatom,jatom,katom,latom,rcd,thbcd,phabcd
            #       MOVE ATOM C TO ORIGIN
            xjl = cart[latom][0] - cart[jatom][0]
            yjl = cart[latom][1] - cart[jatom][1]
            zjl = cart[latom][2] - cart[jatom][2]
            xb = cart[katom][0] - cart[jatom][0]
            yb = cart[katom][1] - cart[jatom][1]
            zb = cart[katom][2] - cart[jatom][2]
            #      ROTATE ABOUT Z-AXIS TO MAKE YB=0, XB POSITIVE. IF XYB IS TOO
            #      SMALL, ROTATE FIRST 90 DEGREES ABOUT THE Y AXIS.
            xyb = math.sqrt(xb * xb + yb * yb)
            k = 1
            if (xyb <= 0.1):
                k = 0
                xpa = zjl
                zpa = -1.0 * xjl
                xjl = xpa
                zjl = zpa
                xpb = zb
                zpb = -1.0 * xb
                xb = xpb
                zb = zpb
                xyb = math.sqrt(xb * xb + yb * yb)
            if (xyb <= epsilon): xyb = epsilon
            costh = xb / xyb
            sinth = math.sqrt(1.0 - costh * costh)
            if (yb <= 0.0): sinth = - sinth
            xpa = xjl * costh + yjl * sinth
            ypa = yjl * costh - xjl * sinth
            #      ROTATE ABOUT THE Y AXIS TO MAKE ZB VANISH.
            rbc = math.sqrt(xb * xb + yb * yb + zb * zb)
            if (rbc <= epsilon): rbc = epsilon
            sinph = zb / rbc
            #      this cosph always positive
            cosph = math.sqrt(1.0 - sinph * sinph)
            xqa = xpa * cosph + zjl * sinph
            zqa = zjl * cosph - xpa * sinph
            #      ROTATE ABOUT X-AXIS TO MAKE ZA=0,YA POSITIVE.
            yza = math.sqrt(ypa * ypa + zqa * zqa)
            if (yza <= epsilon): yza = epsilon
            coskh = ypa / yza
            sinkh = math.sqrt(1.0 - coskh * coskh)
            if (zqa < 0.0): sinkh = -1.0 * sinkh
            #      THE COORDINATES OF A,B, AND C ARE= (XQA,YZA,0),(RBC,0,0) AND
            #      (0,0,0), RESPECTIVELY. THE COORDINATES OF D ARE NOW CALCULATED
            #      IN THE NEW FRAME.
            sina = math.sin(thbcd)
            xd = rcd * math.cos(thbcd)
            yd = rcd * sina * math.cos(phabcd)
            zd = rcd * sina * math.sin(phabcd)
            #      TRANSFORM COORDINATES OF D BACK TO ORIGINAL SYSTEM.
            ypd = yd * coskh - zd * sinkh
            zpd = zd * coskh + yd * sinkh
            xpd = xd * cosph - zpd * sinph
            zqd = zpd * cosph + xd * sinph
            xqd = xpd * costh - ypd * sinth
            yqd = ypd * costh + xpd * sinth
            if (k != 1):
                xrd = -1.0 * zqd
                zrd = xqd
                xqd = xrd
                zqd = zrd
            temp = [xqd + cart[jatom][0], yqd + cart[jatom][1], zqd + cart[jatom][2]]
            cart[iatom] = temp;

    cart = cart[3:len(cart)]
    return cart


####################################
def read_zmat_template(filename):
    f = open(filename, "r")
    zmat = []
    parent = []
    names = []
    line = f.readline()
    while f:  #Get past coments
        if (not re.search('^\*', line)): break
        line = f.readline()
    a = re.search('^\s*\S+\s+(\d+)', line)
    if (a):
        num_atoms = int(a.group(1))
    else:
        raise Exception("Error reading zmat from template \n%s\n%s" % (filename, line))
    for i in range(num_atoms):
        line = f.readline()
        a = re.search('^\s*\d+\s+(\d+)\s+\S+\s+\S+\s+(\S+)\s+\d+\s+([\d\.\-]+)\s+([\d\.\-]+)\s+([\d\-\.]+)', line)
        if (a):
            temp = [float(a.group(3)), float(a.group(4)), float(a.group(5))]
            if (len(zmat) <= 0):
                zmat = [[float(a.group(3)), float(a.group(4)), float(a.group(5))]]
            else:
                zmat.append(temp)
            parent.append(int(a.group(1)) - 1)
            #de moemnt
            names.append(a.group(2).strip('_'))
        else:
            raise Exception("Error reading zmat from template")
    f.close()
    return [parent, zmat, names]


#################################################
def MatchTempMaeAtoms(mae_file, template_file):
    """

      Say which number of atom of th mae_file corresponds to the one on the template_file

      Input:
        mae_file: topology of the ligand
        template_file: topology of the ligand template

      Output:
        mae2temp: Which number of atom of the template file correspond to the same atom on the mae file
        Example-->[0, 1, 2, 3, 4, 6, 8, 10, 12, 13, 14, 15, 16, 5, 7, 9, 11, 17, 18, 19]
        Explanation--> Atom number 5 from the mae correspond the 6th of the template.

        temp2mae: Which number of atom of the mae file correspond to the same atom on the template file
        Example--> [0, 1, 2, 3, 4, 13, 5, 14, 6, 15, 7, 16, 8, 9, 10, 11, 12, 17, 18, 19]
        Explanation--> Atom number 5 of the template file corresponds to atom number 13 of the mae file.

    """
    [parent, zmat, temp_names] = read_zmat_template(template_file)
    mae_names = find_names_in_mae(mae_file)
    if ( len(temp_names) != len(mae_names)):
        raise Exception(
            "MAE and template of different length %d!=%d (check if the input file has more then one residue)\n" % (
            len(mae_names), len(temp_names)));
    mae2temp = []
    temp2mae = []
    num_found = 0
    for i in range(len(temp_names)):
        temp_names[i] = re.sub('_', ' ', temp_names[i])
        mae2temp.append(-100)
        temp2mae.append(-100)
    for i in range(len(mae_names)):
        for j in range(len(temp_names)):
            if (temp_names[j] == mae_names[i] ):
                mae2temp[i] = j
                temp2mae[j] = i
                num_found = num_found + 1
            #    print "Match ","[",mae_names[i],"][",temp_names[mae2temp[i]],"]"
    if ( len(mae2temp) < len(mae_names)):
        raise Exception("Mae file and template do not match\n");
    if (len(mae2temp) != num_found):
        ERROR = "Mae file and template do no match: Check to see if all mae atoms have a pdb atom name\nAtoms not Found:  ";
        for i in range(len(mae_names)):
            if (mae2temp[i] == -100): ERROR += " [";ERROR += mae_names[i];ERROR += "]";
        raise Exception(ERROR);
    #  for i in range(len(mae_names)):
    #     print i,mae_names[i],mae2temp[i],temp_names[mae2temp[i]]
    return [mae2temp, temp2mae]


#################################################
def FindTorsAtom(tors, tors_ring_num, parent):
    '''
        Find atoms with torsion and append them on the zmatrix

        Input:
            tors: atoms with torsion
            tors_ring_num: ring atoms with torsion
            parent: atom connectivity of ordering

        Output:
            out_tors: torsions
            out_tors_ring_num: ring atom torsions
            zmat_atoms: atoms to rotate
    '''

    zmat_atoms = []
    out_tors = []
    out_tors_ring_num = []
    for i in range(len(tors)):
        found = 0
        for iatom in range(len(parent)):
            jatom = parent[iatom];
            katom = parent[jatom]
            if ( (jatom == tors[i][0] and katom == tors[i][1])
                 or (jatom == tors[i][1] and katom == tors[i][0])):
                if (found == 0):
                    zmat_atoms.append(iatom);
                    found = 1
                    out_tors.append(tors[i])
                    if (len(tors_ring_num) > 0):
                        out_tors_ring_num.append(tors_ring_num[i])
                else:
                    if (iatom < zmat_atoms[-1]): 
                      zmat_atoms[-1] = iatom
        if (found == 0): pass  # No dependents, it must be the dependent in a ring
    return [out_tors, out_tors_ring_num, zmat_atoms]


################################
def get_opts(argv):
    opts = [];
    inputs = []
    for i in argv[1:]:
        a = re.search(r'^[-](.*)', i)
        if (a):
            opts.append(a.group(1))
        else:
            inputs.append(i)
    return inputs, opts


#################################################
def get_first_res_name(filename):
    f = open(filename, 'r')
    res_name = ''
    while 1:
        line = f.readline()
        if line == "": break
        a = re.search(r'^[HA][ET][TO][AM][T ][M ].{11}(...).(.)\s*(\d+)', line)
        if (a):
            c = a.group(2)
            if (c == ' '): c = '_'
            res_name = c + ":" + a.group(3)
            res_type = a.group(1)
            break;
    f.close()
    return [res_name, res_type]


#################################################
#def make_plop_control(conf_file,first_pdb,template,res_name,res_type,libname,
#    gridres,zmat_atoms,use_rings):
#  if(libname==""):
#    temp_file_name = "make_rotamer_library_"+res_type+".con"
#    temp_log_file =  "make_rotamer_library_"+res_type+".log"
#    temp_out_file =  "make_rotamer_library_"+res_type+".out"
#  else:
#    temp_file_name = "make_rotamer_library_"+libname+".con"
#    temp_log_file =  "make_rotamer_library_"+libname+".log"
#    temp_out_file =  "make_rotamer_library_"+libname+".out"
#  if(use_rings): gridres='1.0'
#  fout=open(temp_file_name,"w")
#  line='file datadir '+plop_data_dir+'\n';        fout.write(line)
#  line='file logfile '+temp_log_file+'\n\n';      fout.write(line)
#  if(first_pdb != ""):
#     line = 'load pdb '+first_pdb+' het yes &\n';      fout.write(line)
#     line='  template '+res_type+' '+template+' \n\n'; fout.write(line)
#  line='rot read '+res_name+' '+conf_file+' &\n'; fout.write(line)
#  line='  gridres '+gridres+' &\n';               fout.write(line)
#  line='  template '+res_type+' '+template+' &\n';fout.write(line)
#  if(libname!=""):
#    line='  libname '+libname+' &\n';              fout.write(line)
#  if(use_rings):
#    line='  downsample no &\n';                   fout.write(line)
#  for atom in zmat_atoms:
#    line='  atom '+atom+' &\n';                   fout.write(line)
#  line='\n\n';                                    fout.write(line)
#  fout.close()
#  line='echo '+temp_file_name+' | '+plop_executable + ' > ' + temp_out_file
#  if(not debug):os.system(line)
#  return temp_file_name

####################################

def find_resnames_in_mae(filename):
    builder = MaeFileBuilder()
    mae = builder.build(filename)

    #Get residue names list, remove duplicates building a set, and return a list with unique values
    return list(set([atom['s_m_pdb_residue_name'].strip() for atom in mae.atoms()]))


###################################



class MaeFile:
    def __init__(self, atoms):
        self.__atoms = atoms

    def atoms(self):
        return self.__atoms


    def __str__(self):
        return '\n'.join([str(atom) for atom in self.__atoms])


class MaeFileBuilder:
    def build(self, fileName):
        f = open(fileName, 'r')
        atoms = self.__buildAtoms(self.__getAtomsData(f))
        f.close()
        return MaeFile(atoms)


    def __getAtomsData(self, file):

        numberOfAtoms = self.__skipUntilAtomsLine(file)
        atomsProperties = self.__collectAtomsProperties(file)
        atomsData = self.__collectAtomsData(file, numberOfAtoms)

        return {
            "number": numberOfAtoms,
            "propertyNames": atomsProperties,
            "data": atomsData
        }


    def __buildAtoms(self, allData):
        ret = []
        for atomRegister in allData["data"]:
            ret.append(self.__buildAtom(allData["propertyNames"], atomRegister))

        return ret


    def __skipUntilAtomsLine(self, file):
        atomsBlockSeen = False
        numberOfAtoms = -1

        while not atomsBlockSeen:
            line = file.readline()
            if 'm_atom' in line:
                atomsBlockSeen = True
                m = re.search('\[(\d+)\]', line)
                numberOfAtoms = int(m.groups()[0])

        if numberOfAtoms == -1:
            raise Exception('Error in mae file format: number of atoms not indicated as expected')

        return numberOfAtoms


    def __collectAtomsProperties(self, file):
        atomsProperties = []
        endOfAtomsProperties = False
        while not endOfAtomsProperties:
            line = file.readline().lstrip()
            if line[0] == '#':
                continue
            elif ':::' in line:
                endOfAtomsProperties = True
            else:
                atomsProperties.append(line.strip('\n'))

        return atomsProperties


    def __collectAtomsData(self, file, expectedNumberOfAtoms):
        ret = []
        for i in range(1, expectedNumberOfAtoms + 1):
            line = file.readline()
            if ':::' in line:
                raise Exception('Actual number of atoms is not the same as indicated in m_atoms header ' \
                                '(Expected ' + str(expectedNumberOfAtoms) + ' Obtained ' + str(i) + ')')
            ret.append(self.__tokenizeLine(line))

        if not ':::' in file.readline():
            print('WARNING: Actual number of atoms is greater than the indicated in m_atoms header ')

        return ret

    def __tokenizeLine(self, line):
        output = []
        while (len(line) > 0):
            a = re.search(r'^\s*(\S+)(.*)', line)
            if (a):
                b = re.search(r'\"', a.group(1))
                if (b):
                    a = re.search(r'^\s*\"([^\"]*)\"(.*)', line)
                    if (a):
                        output.append(a.group(1))
                        line = a.group(2)
                    else:
                        print("Error in mae formating {}\n".format(line))
                        return -1;
                else:
                    output.append(a.group(1))
                    line = a.group(2)
            a = re.search(r'^\s*$', line)
            if (a): break;
        return output


    def __buildAtom(self, propertyNames, values):
        atom = {'id': values[0]}
        atomId = values[0]

        if len(propertyNames) + 1 != len(values):
            raise Exception('Atom with id ' + str(atomId) + ' does not match the specified properties')

        for i in range(1, len(values)):
            atom[propertyNames[i - 1]] = values[i]

        return atom
############################################################################

def build_template(mae_file, output_template_file):
  """
    Build ligand template from mae file

    Input:
      mae_file: input .mae
      outpt_template_file = output path

    Output:
      output: ligand template
      resname: residue name

    Explanation:

    Preparation of param:
      1- Use the flld_server schrodinger command line server
         to create the 2005OPLS paramaters for the ligand
      2- Parse the previous parameters extracting 
         atomtypes, bond info, non bonded info, torsions,
         dihedrals and imporpers rotations.
      3- From this parameters create the internal coordinates
         and get some extra parameters as the nonbonded
         sgb_radius, gammas, alphas from our proper database
         which is in PARAM_PATH in case it needs to be extended
         by the user.
      4- With the previous parsed parameters create the several
         sections of the template
      5- Write everything to a RESIDUENAME.hetgrp file
  """
  
  #Create ligand params with ffld_sever from schrodinger
  ffld_server_command = os.path.join(os.environ['SCHRODINGER'], 'utilities/ffld_server')
  subprocess.call([ffld_server_command, "-imae", mae_file, "-version",
    OPLS_VERSION, "-print_parameters", "-out_file", OPLS_CONVERSION_FILE])

  #Retrieve all the useful information from that params
  atom_names_param = retrieve_atom_names(OPLS_CONVERSION_FILE)
  search_and_replace(OPLS_CONVERSION_FILE, atom_names_param)
  atom_types, parents, charges, sigmas, epsilons, stretchings, tors, phis, impropers = parse_param(OPLS_CONVERSION_FILE, len(atom_names_param))

  #Connectivity information from Mae
  res_name = find_resnames_in_mae(mae_file)[0] #Ligand must be defined as a whole residue 
  atom_names = find_names_in_mae(mae_file, undersc=True)
  # bonds = [[stretching[0:2] for stretching in stretchings]
  zmat = create_zmatrix(mae_file, parents)
  number_bonds = len(stretchings)
  number_torsions = len(tors)
  number_atoms = len(atom_names)
  number_phis = len(phis)
  number_improper = len(impropers)

  #Missing NBND information in param
  sgb_radius, gammas, alphas = SGB_paramaters(atom_types)

  #Output file
  output_file = output_template_file + '.hetgrp_ffgen'  
  ################################Template Creation########################33
  
  header = ["* LIGAND DATABASE FILE (OPLS2005)",
            "*",
            '{0:>0}  {1:>6} {2:>5} {3:>5} {4:>7} {5:>7}'.format(res_name, number_atoms, number_bonds,
                                                      number_torsions, (number_phis + number_improper), 0)
            ]

  connectivity_section = []
  for i, (atom_name, atom_type, parent, zmat_row) in enumerate(
    zip(atom_names, atom_types, parents, zmat), 1):
    connectivity_line = '{0:>5} {1:>5} {2:>0} {3:>5} {4:>5} {5:>5} {6:>11.5f} {7:>11.5f} {8:>11.5f} '.format(
      i, parent, 'S', atom_type, atom_name,
      UNK_INT, zmat_row[0], zmat_row[1], zmat_row[2])
    connectivity_section.append(connectivity_line)

  NBOND_section = []
  for i, (atom_type, charge, sigma, epsilon, radius, gamma, alpha) in enumerate(
    zip(atom_types, charges, sigmas, epsilons, sgb_radius, gammas, alphas), 1):
        NBOND_section.append('{0:>5} {1:>8.4f} {2:>8.4f} {3:>10.6f} {4:>8.4f} {5:>8.4f} {5:>13.9f} {6:>13.9f}'.format(
          i, float(sigma), float(epsilon), float(charge), float(radius), float(sigma)/2.0, float(gamma), float(alpha)))
        #  NBOND_section.append('{0:>5} {2:>6} {9:>6} {1:>8} {10:>6} {11:>6} {3:>11.8f} {4:>11.8f}'.format(i, charge, sigma,epsilon, vdw))
 
  strech_section = []
  for stretching in stretchings:
        strech_section.append('{0:>5} {1:>5} {2:>9.3f}  {3:>5.3f}'.format(
          stretching[0], stretching[1], float(stretching[2]), float(stretching[3])))
   
  tors_section = []
  for tor in tors:
    tors_section.append('{0:>5} {1:>5} {2:>5} {3:>11.5f} {4:>10.5f}'.format(
      tor[0], tor[1], tor[2], float(tor[3]), float(tor[4])))

  phi_section = []
  phis = descompose_dihedrals(phis)
  phis = negative_torsions_for_pele(phis)
  for phi in phis:
    phi_section.append('{0:>5} {1:>5} {2:>5} {3:>5} {4:>9.5f} {5:>4.1f} {6:>3.1f}'.format(
      phi[0], phi[1], phi[2], phi[3], (float(phi[4])/2.0), phi[5], abs(float(phi[6]))))

  iphi_section = []
  for improper in impropers:
    iphi_section.append('{0:>5} {1:>5} {2:>5} {3:>5} {4:>9.5f} {5:>3.1f} {6:>3.1f}'.format(
      improper[0], improper[1], improper[2], improper[3], float(improper[4])/2.0, -1, 2))

  file_content = []
  file_content.extend( header +
                       connectivity_section +
                       ['NBON'] +
                       NBOND_section +
                       ['BOND'] +
                       strech_section +
                       ['THET'] +
                       tors_section +
                       ['PHI'] +
                       phi_section +
                       ['IPHI'] +
                       iphi_section +
                       ['END'])

  #Write to file
  with open(output_file, 'w') as f:
    f.write('\n'.join(file_content))
  
  #Remove param.dat file
  try:
    os.remove(OPLS_CONVERSION_FILE)
  except OSError:
    print("Error, param.dat not created. Make sure $SCHRODINGER/utilities/ffld_server is up and running in your computer.")

  print("Template {} generated successfully".format(output_file))
  return output_file, res_name, mae_file, output_file, res_name


def SGB_paramaters(atom_types, tried = []):
  radius = []
  gammas = []
  alphas = [] 
  lines = preproces_file_lines(PARAM_PATH)
  for atom_type in atom_types:
    found = False
    for line in lines:
      line = line.split()
      atom_type_file = line[1]
      if(atom_type == atom_type_file):
        radius.append(line[4])
        gammas.append(line[6])
        alphas.append(line[7])
        found = True
    if not found:
      new_params = find_similar_atomtype_params(atom_type, tried=[])
      if(new_params):
        radius.append(new_params[4])
        gammas.append(new_params[6])
        alphas.append(new_params[7])
        warnings.warn("Paramaters of {} used for {}.".format(
          new_params[1], atom_type))
      else:
        raise Exception('ATOMTYPE PARAMATERS for [{}] NOT FOUND. \
Include them in PlpRotTemp/param/sgbnp.param and run \
the program again.'.format(atom_type))
  return(radius, gammas, alphas)



def find_similar_atomtype_params(atom_type, tried):
  """
    Find parameters from a similar atomtype
    according to the table in param/similarity.dat
  """
  new_atom_type=False
  similarity = 0

  #Search for the most similar atom
  lines = preproces_file_lines(SIMILARITY_PATH)
  for line in lines:
    line = line.split()
    if(line[0]==atom_type and line[2]>similarity and line[1] not in tried):
      similarity = line[2]
      new_atom_type = line[1]
    elif(line[1]==atom_type and line[2]>similarity  and line[0] not in tried):
      similarity = line[2]
      new_atom_type = line[0]
  if new_atom_type:
    tried.append(new_atom_type)

    lines = preproces_file_lines(PARAM_PATH)
    for line in lines:
      line = line.split()
      atom_type_file = line[1]
      if(atom_type == atom_type_file):
        return(line)
    new_params = find_similar_atomtype_params(new_atom_type, tried)
    return new_params

  else:
    return []

def create_zmatrix(mae_file, parents):
  """
    Retrieve the internal coordinates (zmatrix)
    from the complex cooridnates & topology:

    Input:
      str1: structure.structure class object
      parents: atom's parent list

    Ouput:
      zmatrix: structure internal coord
  """
  str1 = structure.StructureReader(mae_file).next()
  order = [i for i in range(len(parents))]
  #Gap of 1 between parents list and the one needed for creating zmat
  parents = [int(parent) - 1 for parent in parents] 
  coordinates = [atom._getXYZ() for atom in str1.atom]
  zmat = xyz2int(coordinates, order, parents)
  return zmat


  
def retrieve_atom_names(OPLS_CONVERSION_FILE):
  """
    parse the param.data and retrieve
    the atom names of the molecule
  """
  atom_names = []
  lines = preproces_file_lines(OPLS_CONVERSION_FILE)
  start_found_NBND = False
  for line in lines:
      if not line.startswith("----------"):
        if(line.startswith("atom type vdw")):
          start_found_NBND = True
        elif(start_found_NBND):
          try:
            line = line.split()
            atom_names.append(line[0])
          except IndexError:
            return atom_names


def preproces_file_lines(file):
  with open(file, "r") as f:
    lines = f.readlines()
    for i, line in enumerate(lines):
      line = re.sub(' +',' ',line)
      line = line.strip('\n').strip()
      lines[i] = line
    return lines



def search_and_replace(file, to_search):
  """
    Search and replace atom_names for numbers
  """
  to_replace = range(1, len(to_search)+1)

  f = open(file,'r')
  filedata = f.read()
  f.close()

  for (item_to_search, item_to_replace) in zip(to_search, to_replace):
      filedata = filedata.replace(' ' + item_to_search.strip('_')+' ', ' '+str(item_to_replace)+' ') #atom types are like _O1_ (strip)

  f = open(file,'w')
  f.write(filedata)
  f.close()

def parse_param(param_file, natoms):
  """
    Parse the OPLS conversion param file
    and get the atomtypes.

    Go through the param.dat through several stages
    while jumping al the ----- lines:


    1. Look for the BSCI section for each atom
       and its parent appending both like
       bond [A,B] and parent [B].
    2. Set end_connectivity_found=True to start
       parsing the NBND section
    
    In the NBND section:
    2. Look for the atom type vd2w line
    3. Parse each line of the non bonding section
    4. set NBND_finished= True to start parsing the bonded section
    

    In the BOND section
    5. For each subsection ([stretchings, bendings, proper_tors, improper_tors]):
       1. Search for the respective title indicating the beggining of the section
       (["Stretch", "Bending", "proper Torsion", "improper Torsion"])
       2. Then parse the section obtaining the values from the columns indicated in columns to take
       3.when all sections are parsed return all values
  """

  atom_types = []
  parents = [0] * natoms
  charges = []
  sigmas = []
  epsilons = []
  vdws = []
  stretchings = []
  bendings = []
  proper_tors = []
  improper_tors = []
  lists = [stretchings, bendings, proper_tors, improper_tors]
  keywords = ["Stretch", "Bending", "proper Torsion", "improper Torsion"]
  columns_to_take = [[0,1, 2, 3], [0, 1, 2, 3,4], [0, 1 ,2 , 3, 4, 5, 6, 7], [0, 1, 2, 3, 4]]
  start_connectivity_found = False
  end_connectivity_found = False
  start_found_NBND = False
  NBND_finished = False
  start_found_BND = False


  with open(param_file) as f:
    lines = f.readlines()
    for i, line in enumerate(lines):
      #prepared lines
      line = re.sub(' +',' ',line)
      line = line.strip('\n').strip()

      if not end_connectivity_found:
        if(line.startswith("BCI's")):
          start_connectivity_found = True
        elif(start_connectivity_found):
          try:
            line = line.split()
            #parents[int(line[1])-1]-->Atom names start at 1 but list at 0
            parents[int(line[1])-1] = int(line[0]) 
          except IndexError:
            end_connectivity_found = True


      #NBND section
      elif (line.startswith("----------") is False) and (NBND_finished is False) and (end_connectivity_found is True):
        if(line.startswith("atom type vdw")):
          start_found_NBND = True
        elif(start_found_NBND):
          try:
            line = line.split()
            atom_types.append(line[3])
            charges.append(line[4])
            sigmas.append(line[5])
            epsilons.append(line[6])
          except IndexError:
            NBND_finished=True

      elif(NBND_finished):
        for (keyword, indexes, List) in zip(keywords, columns_to_take, lists):
          while(line.startswith(keyword) is False):
            line, i = move_line_forward(lines, i)
          line, i = move_line_forward(lines, i)
          while(line):
                line = re.sub(' +',' ',line)
                line = line.strip('\n')
                line = line.split()
                values = [line[index] for index in indexes]
                List.append(values)
                line, i = move_line_forward(lines, i)
        return atom_types, parents, charges, sigmas, epsilons, stretchings, bendings, proper_tors, improper_tors

def move_line_forward(lines, i):
  i+=1
  line = lines[i].strip()
  return line, i

def descompose_dihedrals(phis):
  """
    For each dihedral line as:
    atom1 atom2 atom3 atom4   V1    V2   V3    V4
      1     2     3      4   0.00 1.000 5.000 3.000

    Separate all components as next:
    atom1 atom2 atom3 atom4   value  which VX???
     1     2     3      4     1.000      2
     1     2     3      4     5.000      3
  """
  new_phis = []
  new_phi = []
  for phi in phis:
    atoms = phi[0:4]
    if(phi[4:8] == ['0.000','0.000','0.000','0.000']):
      new_phis.append([atoms[0], atoms[1], atoms[2], atoms[3], phi[4], 1, 1])
      continue
    for index, component in enumerate(phi[4:8], 1):
      if(component != '0.000'):
        new_phis.append([atoms[0], atoms[1], atoms[2], atoms[3], component, 1, index])
  return new_phis

def negative_torsions_for_pele(phis):
  for i, phi in enumerate(phis):
    atom1 = phi[0]
    atom4 = phi[3]
    for phi_to_compare in phis[0:i]:
      
      if(atom1 == phi_to_compare[0] and atom4 == phi_to_compare[3] and phi_to_compare[2]<0):
        break
      elif(phi[0:4] == phi_to_compare[0:4]):
        pass
      elif(phi_to_compare[2]<0):
        pass
      elif(atom1 in phi_to_compare and atom4 in phi_to_compare and phi!=phi_to_compare):
        phi[2] = -int(phi[2])
    if(int(phi[6]) in [1,3]):
      phi[5] = 1.0
    else:
      phi[5] = -1.0
  return(phis)





##################################################
def make_lib_from_mae(lib_name, lib_type, conf_file, tors, names, parent, ordering, mae2temp, temp2mae, gridres,
                      lib_gridres):
    allow_downsample = 1
    if (float(lib_gridres) < 0):
        lib_gridres = float(lib_gridres) * -1
        allow_downsample = -1
    # Get The ZMatrix atom for each tors
    [tors, holder, zmat_atoms] = FindTorsAtom(tors, [], parent)
    # Convert back to atom numbering in mae file
    zmat_atoms_mae = []
    for i in range(len(zmat_atoms)):
        zmat_atoms_mae.append(temp2mae[ordering[zmat_atoms[i]]])

    # for i in range(len(zmat_atoms)):
    #    print "%d %s %d %d"%(zmat_atoms[i],names[zmat_atoms[i]],zmat_atoms_mae[i],ordering[zmat_atoms[i]]);

    ordering_to_mae = []
    for i in range(len(ordering)):
        ordering_to_mae.append(temp2mae[ordering[i]])

    # Read in File
    ###############################CHANGE SCHRODINGER########################
    st1 = structure.StructureReader(conf_file).next()
    tors_values = []
    for st in structure.StructureReader(conf_file):
        # Read in the xyz coordinates into cart
        cart = [atom._getXYZ() for atom in str1.atom]

        #    print 'cartesian(xyz)'
        #    for i in range(len(cart)):
        #       print 'X',cart[i][0],cart[i][1],cart[i][2]
        zmat = xyz2int(cart, ordering_to_mae, parent)
        this_tors = []
        #    print 'zmat'
        #    for i in range(len(zmat)):
        #       print 'X',zmat[i][0],zmat[i][1],zmat[i][2]
        for atom in zmat_atoms:
            this_tors.append(zmat[atom][2])
        tors_values.append(this_tors)
    ###############################CHANGE SCHRODINGER########################
    # Write To Library File
    # for atom in zmat_atoms:
    #    print "%d %s\n"%(atom,names[atom]); 
    # for tv in tors_values:
    #    print tv,"\n";
    # Change to Grid Resolution
    grid_entries = []
    #print "GRIDRES %f"%float(gridres)
    if (float(gridres) != float(lib_gridres)):
        unique_tors_values = []
        for i in range(len(tors_values)):
            min_diff = float(gridres) * 1000
            for j in range(len(unique_tors_values)):
                diff = 0.0
                for k in range(len(tors_values[i])):
                    temp = [abs(tors_values[i][k] - unique_tors_values[j][k]), \
                            abs(tors_values[i][k] - unique_tors_values[j][k] + 360), \
                            abs(tors_values[i][k] - unique_tors_values[j][k] - 360)]
                    if (min(temp) > diff): diff = min(temp)
                if (diff < min_diff): min_diff = diff
            if (min_diff > float(gridres)): unique_tors_values.append(tors_values[i])
        print("Reduced number of conformers from {0} to {1} at resolution  {2}".format(len(tors_values),len(unique_tors_values),str(gridres)))
        tors_values = unique_tors_values

    for tv in tors_values:
        temp = [];
        for element in tv:
            element = round(element / float(lib_gridres))
            temp.append(element)
        grid_entries.append(temp)
    #for ge in grid_entries:
    #   print ge;
    # Sort the entries
    for i in range(len(grid_entries)):
        for j in range(len(grid_entries)):
            if (grid_entries[i] < grid_entries[j]):
                temp = grid_entries[i]
                grid_entries[i] = grid_entries[j]
                grid_entries[j] = temp
                #print "Sorted"
    #for ge in grid_entries:
    #   print ge;
    # Eliminate Duplicates
    unique_grid_entries = [];
    for i in range(len(grid_entries)):
        if (i == 0 or grid_entries[i - 1] != grid_entries[i]):
            unique_grid_entries.append(grid_entries[i])
    grid_entries = unique_grid_entries

    #print "Unique Sorted"
    #for ge in unique_grid_entries:
    #   print ge;

    #Actually Write The Library File
    fp = open(lib_name + "." + lib_type, "w")
    line = "* %6s %2d %7d %4.1f 0\n" % (
    lib_name, len(zmat_atoms), len(grid_entries), float(lib_gridres) * allow_downsample)
    fp.write(line);
    for atom in zmat_atoms:
        fp.write("%4s\n" % names[atom])
    for ge in grid_entries:
        line = "";
        for element in ge:
            line = line + " %5d" % element;
        line = line + "\n";
        fp.write(line)
    fp.close


##################################################
def find_yesno(input):
    a = re.search('[yY][eE][sS]', input);
    if (a): return 1
    a = re.search('[nN][on]', input);
    if (a): return 0
    a = re.search('[Oo][Nn][Ll][Yy]', input);
    if (a): return 2
    a = re.search('[yY]', input);
    if (a): return 1
    a = re.search('[nN]', input);
    if (a): return 0
    a = re.search('[Oo]', input);
    if (a): return 2
    print("ERROR PROCESSING {} as yes/no".format(input))
    return -1


##################################################
def make_libraries(resname, conf_file, root, names, zmat_atoms, group, use_rings, use_mult_lib, output_template_file,
                   gridres, debug, R_group_root_atom_name='None'):
    #Convert conf file to pdbs if necessary
    is_pdb = re.search(r'^(.*)\.[PpEe][DdNn][BbTt]', conf_file)
    is_mae = re.search(r'^(.*)\.[Mm][Aa][Ee]', conf_file)
    if (not is_mae):
        print("Conformation file must be an mae file\n")
        sys.exit(-1)

    assign_filename = resname.upper() + ".rot.assign"
    line = "rot assign res " + resname.upper() + " &\n";
    # print 'template file name: ',output_template_file
    #  if(use_mult_lib==0):
    #    zmat_names=[]
    #    for i in zmat_atoms:
    #      zmat_names.append(names[i])
    #    line='';
    #    for i in zmat_names: line = line + i.rjust(5)
    #    libname=res_type +"__1"
    #    print 'Building Rotamer Library: ',libname+'.side'
    #    print "Zmat Atoms for rotatable bonds:", line
    #    plop_con = make_plop_control(conf_file,pdb_root,output_template_file,
    #      res_name,res_type,libname,gridres,zmat_names,use_rings)
    #    line = "rot assign res "+res_type+" sidelib "+libname+" default\n"; 
    #    print "\nUse this line to utilize the library(",assign_filename,"):\n",line;
    #    fp=open(assign_filename,"w")
    #    fp.write(line);fp.close();
    #  else:

    ######################CHANGE SCHRODINGER#########################
    if use_rings:
        lib_gridres = -1.0  # If there are flexible rings present we set the library to be at res 1 and
    # Disalow down sampling
    else:
        lib_gridres = gridres
    ######################CHANGE SCHRODINGER#########################

    for grp in range(max(group) + 1):
        grp_tors = []
        if (grp + 1 < 10):
            libname = "%3s__%1d" % (resname.upper(), grp + 1)
        else:
            libname = "%3s_%2d" % (resname.upper(), grp + 1)

        still_mainchain = 0
        if (R_group_root_atom_name != 'None'):
            still_mainchain = 1
            R_group_root_atom_order_num = convert_name_to_num([R_group_root_atom_name], names)

        for t in tors:
            torsion_atom = t[1]
            myname = convert_num_to_name([int(torsion_atom)], names)
            myname_formatted = re.sub(' ', '_', str(myname[0]))
            #If we've already gone up to or past the R group, start counting it as the side-chain
            if (R_group_root_atom_name != 'None'):
                if ((myname_formatted == R_group_root_atom_name) or (
                    torsion_atom >= R_group_root_atom_order_num[0])): still_mainchain = 0
            if ( group[t[0]] == grp or group[t[1]] == grp and still_mainchain == 0):
                grp_tors.append(t)
        if (grp != 0 ): line += "  newgrp &\n";
        line += "  sidelib " + libname + " default &\n"
        make_lib_from_mae(libname, "side", conf_file, grp_tors, names, parent, old_atom_num, mae2temp, temp2mae,
                          gridres, lib_gridres)
    fp = open(assign_filename, "w")
    fp.write(line);
    fp.close();
    return assign_filename


################################################################################
def convert_gridres(gridres):
    # Write the information for libaries
    # Get the correct library
    if (eval(gridres) == 5.0):
        lib_name = "FREE_5"
    elif (eval(gridres) == 10.0):
        lib_name = "FREE10"
    elif (eval(gridres) == 15.0):
        lib_name = "FREE15"
    elif (eval(gridres) == 20.0):
        lib_name = "FREE20"
    elif (eval(gridres) == 30.0):
        lib_name = "FREE30"
    elif (eval(gridres) == 40.0):
        lib_name = "FREE40"
    elif (eval(gridres) == 45.0):
        lib_name = "FREE45"
    elif (eval(gridres) == 60.0):
        lib_name = "FREE60"
    elif (eval(gridres) == 90.0):
        lib_name = "FREE90"
    elif (eval(gridres) == 180.0):
        lib_name = "FRE180"
    else:
        raise Exception('Incorrect grid resolution (5,10,15,20,30,40,45,50,60,90,180)')
    return lib_name


################################################################################
def check_oh_tors(mae_file, tors, names):
    # see if each tors is for an OH (or O-D etc) by checking to see if the 
    # larger atom number (farther from core) is bound to only one atom
    # which is in turn bound to no atoms
    natoms = len(names)
    st1 = structure.StructureReader(mae_file).next()
    is_oh = []
    for tor in tors:
        atom_name0 = names[tor[0]].replace("_", " ");
        atom_name1 = names[tor[1]].replace("_", " ");
        mae_atom1 = -1
        for iatom in range(natoms):
            if (atom_name1.strip() == (st1.atom[iatom + 1].property['s_m_pdb_atom_name']).strip()):
                mae_atom1 = iatom + 1
        if (mae_atom1 < 0):
            raise Exception('Error in check OH Tors')
        bound_atoms = []
        for iatom in range(natoms):
            if (atom_name0 != st1.atom[iatom + 1].property['s_m_pdb_atom_name'] and
                        st1.getBond(mae_atom1, iatom + 1) != None):
                bound_atoms.append(iatom + 1)
            #    print "ATOM ",st1.atom[mae_atom1].property['s_m_pdb_atom_name'],len(bound_atoms)
            #    for iatom in bound_atoms:
            #       print "Bound ",st1.atom[iatom].property['s_m_pdb_atom_name']

        if (len(bound_atoms) == 1):
            any_bound = 0;
            first_bound_atom = bound_atoms[0]
            for iatom in range(natoms):
                if (atom_name1 != st1.atom[iatom + 1].property['s_m_pdb_atom_name'] and
                            st1.getBond(first_bound_atom, iatom + 1) != None):
                    #            print "BOUND TO ",st1.atom[iatom+1].property['s_m_pdb_atom_name']
                    any_bound = 1
            if (any_bound == 0):
                is_oh.append(1)
            else:
                is_oh.append(0)
        else:
            is_oh.append(0)
    return is_oh


################################################################################

def find_build_lib(resname, mae_file, root, tors, names, group, gridres, gridres_oh, use_rings, back_lib, tors_ring_num,
                   ring_libs, debug):
    #Make a filename with this data
    assign_filename = resname.upper() + ".rot.assign";
    f = open(assign_filename, "w")
    lib_name_nom = convert_gridres(gridres)
    lib_name_oh = convert_gridres(gridres_oh)
    is_oh = check_oh_tors(mae_file, tors, names)


    #  if(use_rings):
    #    print "Libraries must be corrected for rings"
    written_ring = []
    if (tors_ring_num != []):
        for i in range(max(tors_ring_num) + 1): written_ring.append(0)
    f.write("rot assign res " + resname.upper() + ' &\n')
    if (back_lib != ""):
        f.write(" backlib " + back_lib + " &\n");
    for grp in range(max(group) + 1):
        for i in range(len(tors)):
            if ( group[tors[i][0]] == grp or group[tors[i][1]] == grp):
                if (tors_ring_num[i] == 0):
                    if (is_oh[i] == 1):
                        lib_name = lib_name_oh
                    else:
                        lib_name = lib_name_nom
                    f.write("   sidelib " + lib_name + " " + names[tors[i][0]] + ' ' + names[tors[i][1]] + " &\n")
                else:
                    ring_num = tors_ring_num[i]
                    if (written_ring[ring_num] == 0):
                        f.write("   sidelib " + ring_libs[ring_num - 1] + "  default &\n")
                        written_ring[ring_num] = 1
        if (grp != max(group)):
            f.write("     newgrp &\n")
    f.close();
    return assign_filename
    


#################################################

def get_root_path(mae_file):
  #File Must By Copied to Local Directory First

  if(os.path.dirname(mae_file)):
    file_name = os.path.basename(mae_file)
    root = os.path.splitext(file_name)[0]
  else:
    root = os.path.splitext(mae_file)[0]
  return root

#TetherRotBonds Class: Chris McClendon, Jacobson Group
#originally, the unnatural amino acid stuff was going to be a separate class that used functions and subroutines from PlopRotTemp.py. It turned out to make more sense to add in many extra functions alongside the rest of the functions and make special cases in the main code for unnatural amino acids, easily retaining backwards-compatibility. Perhaps this should be converted stylisticly to parallel the existing code.

class TetherRotBonds:
    maefile = ""
    logfile = ""
    tors = []
    ring_tors = []
    chain = ''
    resno = ''
    template = ''
    num_atoms = 0
    num_bonds = 0
    num_tetha = 0
    names = []
    bonds = []

    def macromodel_torsions(self):
        # the following parameters are from Ken Borelli's PlopRotTemp.py
        # Create ComUtil instance , define potential energy surface: solution phase, OPLSAA
        # Serial mode enabled so each structure is used to seed a unique search
        mcu = mu.ComUtil(ffld='opls2005', serial=True, solv=True, nant=True, demx=True)
        mxu = mu.CluUtil()

        # There are two debug switch associated with AUTO
        # Debug output appears in jobname.log
        mcu.SOLV[2] = 1  # water
        mcu.DEBG[1] = 520
        mcu.DEBG[2] = 521
        mcu.MCMM[1] = 1  # Take 1,000 steps
        mcu.MCMM[2] = 1  # Store up to 1000 values
        mcu.MINI[3] = 0  # Don't minimize
        mcu.DEMX[5] = 100
        mcu.DEMX[6] = 500
        a = re.match(r'(.*)\.mae', self.maefile)
        root = a.group(1)
        conf_root = root + "_conf"
        com_file = mcu.mcmm(self.maefile, conf_root + '.com')
        logfile = conf_root + '.log'
        cmd = mcu.getLaunchCommand(com_file)
        job = jc.launch_job(cmd)
        job.wait()
        return logfile

    def __init__(self, maefile, chain, resno):
        self.maefile = maefile
        self.logfile = self.macromodel_torsions()
        self.tors = find_tors_in_log(self.logfile)  #uses (atom number - 1)
        [self.ring_tors, ring_num] = find_tors_in_rings(self.tors, self.maefile)
        #remove rings for time being
        self.tors = remove_tors(self.tors, self.ring_tors)
        self.chain = chain
        self.resno = resno
        self.names = find_names_in_mae(self.maefile)
        self.bonds = find_bonds_in_mae(self.maefile)

    def __init__(self, maefile, chain, resno, logfile):
        self.maefile = maefile
        self.logfile = logfile
        self.tors = find_tors_in_log(self.logfile)  #uses (atom number - 1)
        [self.ring_tors, ring_num] = find_tors_in_rings(self.tors, self.maefile)
        #remove rings for time being
        self.tors = remove_tors(self.tors, self.ring_tors)
        self.chain = chain
        self.resno = resno
        self.names = find_names_in_mae(self.maefile)
        self.bonds = find_bonds_in_mae(self.maefile)

    def __init__(self, maefile, chain, resno, logfile, mytors):
        self.maefile = maefile
        self.logfile = logfile
        self.tors = mytors
        self.chain = chain
        self.resno = resno
        self.names = find_names_in_mae(self.maefile)
        self.bonds = find_bonds_in_mae(self.maefile)

    def ismethyl(self, atom):
        atom_is_carbon = re.search('\s*C', self.names[atom])
        group_is_methyl = 0
        if (atom_is_carbon):
            connected = {}
            connected = Buildup_Connected(self.bonds, len(self.names), connected)
            num_hydg = 0
            for i in connected:
                is_hydg = re.search('\s*[0-9]*H', self.names[i])
                if (is_hydg): num_hydg = num_hydg + 1
            if (num_hydg == 3): group_is_methyl = 1
        return group_is_methyl

    def output_rotbonds(self, R_group_root_atom_name='None'):
        fout = open(self.maefile + "_torsions.txt", "w")
        print("\n")
        print("Add the following to your PLOP control or input file")
        print("or add the line from " + str(self.maefile + "_torsions.txt"))
        print("\n")
        numtors = 0
        Calpha_atom = FindCAlpha(self.maefile, raise_exception=False)
        C_atom = FindC(self.maefile, raise_exception=False)

        #Keep as mainchain the z-matrix atoms until R_group_root_atom_name,
        # the root atom of the group to be sampled, is reached
        still_mainchain = 0
        #if(R_group_root_atom_name > -1): still_mainchain=1
        #print "R group root atom:" + str(R_group_root_atom_name) + "\n"

        if (R_group_root_atom_name != 'None'):
            still_mainchain = 1
            R_group_root_atom_order_num = convert_name_to_num([R_group_root_atom_name], self.names)

        for mytors in self.tors:
            torsion_atom = mytors[1]
            myname = convert_num_to_name([int(torsion_atom)], self.names)
            myname_formatted = re.sub(' ', '_', str(myname[0]))
            #print "Atom: "+str(torsion_atom)+" Name of atom "+str(myname_formatted)+" root atom name "+str(R_group_root_atom_name)+" root atom num "+str(R_group_root_atom_order_num[0])+"\n"
            if (R_group_root_atom_name != 'None'):
                if ((myname_formatted == R_group_root_atom_name) or torsion_atom >= R_group_root_atom_order_num[
                    0]): still_mainchain = 0

            if ((not (self.ismethyl(torsion_atom))) and (torsion_atom != Calpha_atom) and (torsion_atom != C_atom) and (
                still_mainchain == 0)):
                numtors = numtors + 1
        print("  nonstandard " + self.chain + ":" + str(self.resno) + " " + str(numtors) + " &")
        fout.write("  nonstandard " + self.chain + ":" + str(self.resno) + " " + str(numtors) + " ")

        still_mainchain = 0
        if (R_group_root_atom_name != 'None'): still_mainchain = 1

        for mytors in self.tors:
            torsion_atom = mytors[1]
            myname = convert_num_to_name([int(torsion_atom)], self.names)
            myname_formatted = re.sub(' ', '_', str(myname[0]))
            if (R_group_root_atom_name != 'None'):
                if ((myname_formatted == R_group_root_atom_name) or torsion_atom >= R_group_root_atom_order_num[
                    0]): still_mainchain = 0
            if ((not (self.ismethyl(torsion_atom))) and (torsion_atom != Calpha_atom) and (torsion_atom != C_atom) and (
                still_mainchain == 0)):
                #myname = convert_num_to_name([int(torsion_atom)],self.names)
                #myname_formatted = re.sub(' ','_',str(myname[0]))
                print("   " + self.chain + ":" + str(self.resno) + ":" + str(myname_formatted) + " &")
                fout.write("   " + self.chain + ":" + str(self.resno) + ":" + str(myname_formatted) + " ")
                #print("\&\n")
        print("\n")
        fout.write("&")
        fout.close()
        print("Use the previous syntax under tether pred, use the following syntax for other modules.")
        return


################################################################################

def build_ring_libs(mae_min_file, root, resname, tors, tors_ring_num, \
                    names, rank, parent, old_atom_num, mae2temp, gridres, files2clean, debug):
    # To do: Prune the ligand so a smaller region has to be sampled.  Currently the 
    # entire ligand is sampled to get the sampling of a single ring.
    nsamp = 1000;
    nrot = 1000
    ring_lib_names = []
    for ring_num in range(1, max(tors_ring_num) + 1):
        libname = resname.upper() + "_R%1d" % ring_num
        conf_root = root + ".ring." + str(ring_num)
        ring_tors = []
        for i in range(len(tors)):
            # Remove All Torsions that are in frozen core areas 
            if (tors_ring_num[i] == ring_num and (rank[tors[i][0]] > 0 or rank[tors[i][1]] > 0)):
                ring_tors.append(tors[i])
        print("Sampling ring ", ring_num, " using macromodel CGEN")
        print("Torsions in ring {}".format(ring_num))
        for t in ring_tors: 
          print ("%s  -  %s" % (names[t[0]], names[t[1]]))
        mcu_conf = mu.ComUtil(ffld='opls2005', serial=True, solv=True, nant=False, demx=True)
        mcu_conf.SOLV[2] = 1  # water
        #     mcu_conf.MCMM[1] = nsamp # Take X steps
        #     mcu_conf.MCMM[2] = nrot # Store up to Y values
        mcu_conf.CGEN[1] = nsamp
        mcu_conf.CGEN[2] = nrot
        mcu_conf.MINI[1] = 1  # PRCG
        mcu_conf.MINI[3] = 50  # iterations of minimze
        #Read in the mae file
        st1 = structure.StructureReader(mae_min_file).next()
        com_file = mcu_conf.mcmm(mae_min_file, conf_root + '.com')
        conf_file = conf_root + '-out.mae'
        #Put constraints on all bonds 
        com_fh = open(com_file, "a")
        ibond = 0
        for i in range(1, len(st1.atom) + 1):
            for j in range(1, len(st1.atom) + 1):
                if (j > i and st1.getBond(st1.atom[i], st1.atom[j]) != None):
                    dist = st1.measure(st1.atom[i], st1.atom[j])
                    mcu_conf.FXDI[8 * ibond + 1] = i
                    mcu_conf.FXDI[8 * ibond + 2] = j
                    mcu_conf.FXDI[8 * ibond + 5] = 99999.99
                    mcu_conf.FXDI[8 * ibond + 6] = dist
                    ibond = ibond + 1
                    opcd_data = mcu_conf.getOpcdArgs('FXDI')
                    com_fh.write(opcd_data)

        #Put constraints on all angles 
        iang = 0
        for i in range(1, len(st1.atom) + 1):
            for j in range(1, i):
                for k in range(i + 1, len(st1.atom) + 1):
                    if (st1.getBond(st1.atom[i], st1.atom[j]) != None and \
                                    st1.getBond(st1.atom[i], st1.atom[k]) != None):
                        angle = st1.measure(st1.atom[j], st1.atom[i], st1.atom[k])
                        mcu_conf.FXBA[8 * iang + 1] = j
                        mcu_conf.FXBA[8 * iang + 2] = i
                        mcu_conf.FXBA[8 * iang + 3] = k
                        mcu_conf.FXBA[8 * iang + 5] = 99999.99
                        mcu_conf.FXBA[8 * iang + 6] = angle
                        iang = iang + 1
                        opcd_data = mcu_conf.getOpcdArgs('FXBA')
                        com_fh.write(opcd_data)
        com_fh.close()

        # Run the search
        if (not debug):
            cmd = mcu_conf.getLaunchCommand(com_file)
            job = jc.launch_job(cmd)
            job.wait()
            files2clean.append(conf_root + '.com')
            files2clean.append(conf_root + '.log')
            files2clean.append(conf_root + '-out.mae')
            files2clean.append(conf_root + '-out.ouL')

        # Build the library for this ring
        make_lib_from_mae(libname, "side", conf_file, ring_tors, names, parent, old_atom_num, mae2temp, temp2mae,
                          gridres, -1.0)
        ring_lib_names.append(libname)

    return ring_lib_names

def printGlobalVariables():
    d = {}
    d["mae_file"] = mae_file
    d["template_file"] = template_file
    d["debug"] = debug
    d["conf_file"] = conf_file
    d["output_template_file"] = output_template_file
    d["gridres"] = gridres
    d["nsamp"] = nsamp
    d["nrot"] = nrot
    d["Ecut"] = Ecut
    d["use_rings"] = use_rings
    d["do_init_min"] = do_init_min
    d["user_core_atom"] = user_core_atom
    d["max_dist_eq"] = max_dist_eq
    d["user_tors"] = user_tors
    d["back_tors"] = back_tors
    d["back_algorithm"] = back_algorithm
    d["back_conf_file"] = back_conf_file
    d["hetgrp_opt"] = hetgrp_opt
    d["use_mae_charges"] = use_mae_charges
    d["OPLS"] = OPLS
    d["max_tors"] = max_tors
    d["user_fixed_bonds"] = user_fixed_bonds
    d["files2clean"] = files2clean
    d["use_mult_lib"] = use_mult_lib
    d["run_conf"] = run_conf
    d["algorithm"] = algorithm
    d["clean"] = clean
    d["gridres_oh"] = gridres_oh
    d["unnat_res"] = unnat_res
    d["resno"] = resno
    d["chain"] = chain
    d["grow"] = grow
    d["tree"] = tree
    d["R_group_root_atom_name"] = R_group_root_atom_name

    print(d)
