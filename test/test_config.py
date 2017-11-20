import sys
import os
import re
import subprocess

NAMES = ['sigma', 'epsilon', 'charges', 'SGBr', 'vdwr', 'alpha', 'gamma']
DEFAULT_RES_NAME = "LIG"
REPOS_PATH = os.path.abspath(os.path.join(__file__ ,"../.."))
TEST_PATH = os.path.join(os.path.dirname(__file__), 'data')
MAIN_PATH =  os.path.join(REPOS_PATH, 'PlpRotTemp/main.py')
OLD_MAIN_PATH = '/home/dani/repos/presentation/PlpRotTemp/main_16_9_17.py'
TRESHOLD= 0.15


try:
    PYTHON_PATH = os.path.join(os.environ['SCHRODINGER'] + "/utilities/python")
except KeyError:
    print("Set SCHRODINGER environment variable path")

class PlopTest:

  def __init__(self, input_file):
    self.input_file = input_file
    self.old_names = None
    self.new_names = None

  def run_versions(self):

    res_name = "LIG"

    try:
      os.remove("LIG".lower())
      os.remove("LIG".upper())
    except OSError:
      pass

    subprocess.call([PYTHON_PATH, MAIN_PATH, os.path.join(TEST_PATH, self.input_file)])
    subprocess.call([PYTHON_PATH, OLD_MAIN_PATH, os.path.join(TEST_PATH, self.input_file)])


    new_sections = self.parse_template(os.path.join(REPOS_PATH,res_name.upper()))
    old_sections = self.parse_template(os.path.join(REPOS_PATH,res_name.lower()))

    self.old_names = old_sections[0][3]
    self.new_names = new_sections[0][3]

    return old_sections, new_sections

  def parse_template(self, file_path):

    res_name = find_resnames_in_mae(os.path.join(TEST_PATH, self.input_file))

    parents, sidechain, atom_types, names = ([] for i in range(4))
    sigma, epsilon, charges, SGBr, vdwr, alpha, gamma = ([] for i in range(7))
    atom1b, atom2b, springK, eq_dist = ([] for i in range(4))
    atom1t, atom2t, atom3t, angles, force = ([] for i in range(5))
    atom1p, atom2p, atom3p, atom4p, valuep, coefp, componentp = ([] for i in range(7))
    atom1i, atom2i, atom3i, atom4i, valuei, coefi, componenti = ([] for i in range(7))

    ITEMS = [
            [parents, sidechain, atom_types, names], 
            [sigma, epsilon, charges, SGBr, vdwr, alpha, gamma],
            [atom1b, atom2b, springK, eq_dist],
            [atom1t, atom2t, atom3t, angles, force],
            [atom1p, atom2p, atom3p, atom4p, valuep, coefp, componentp],
            [atom1i, atom2i, atom3i, atom4i, valuei, coefi, componenti]
            ]

    INDEXES = [
              [1,2,3,4],
              [1,2,3,4,5,6,7],
              [0,1,2,3],
              [0,1,2,3,4],
              [0,1,2,3,4,5,6],
              [0,1,2,3,4,5,6]
              ]
    
    KEYWORDS = [res_name, 'NBON', 'BOND', 'THET', 'PHI', 'IPHI']
    
    starts = [False, False, False, False, False, False]
    ends = [False, False, False, False, False, False]
    print(file_path)
    with open(file_path, 'r') as f:
      line = f.readline() #Skip Commentaris
      line = f.readline() #Skip Commentaris
      line = f.readline() #Skip Commentaris
      for i, (keyword, item, index)  in enumerate(zip(KEYWORDS, ITEMS, INDEXES)):
          line = re.sub(' +',' ',line)
          line = line.strip('\n').strip()
          while not ends[i]:
              if line.startswith(keyword):
                starts[i] = True
              elif starts[i]:
                line = line.split()
                if len(line)==1:
                  ends[i] = True
                  line = line[0]
                  continue
                else:
                  for k, j in zip(item, index):
                      k.append(line[j])
                
              line = f.readline()
    return ITEMS

  
  def compare_NBND_sections(self, old_sections, new_sections):

      num_atoms = len(self.old_names)
    
      NBND_old_sect = old_sections[1]
      NBND_new_sect = new_sections[1]

      NBND_new_ordered = []
      for i in range(num_atoms):
        old_name = self.old_names[i].strip('_')
        for j in range(num_atoms):
          if(old_name == self.new_names[j]):
            NBND_new_ordered.append([param[j] for param in NBND_new_sect])

      translated_old_NBND = zip(*NBND_old_sect)

      return translated_old_NBND, NBND_new_ordered

  def compare_BND_sections(self, old_sections, new_sections, section):


    SECTIONS = { 
                 'BND' : ['BND', 2, 2],
                 'TOR' : ['TOR', 3, 3],
                 'PHI' : ['PHI', 4, 4],
                 'IPHI' : ['IPHI', 5, 4]
              }

    section = SECTIONS[section]
    section_name = section[0]
    section_index = section[1]
    section_elements = section[2]

    BND_old_sect = old_sections[section_index]
    BND_new_sect = new_sections[section_index]

    if(len(BND_old_sect)!=BND_new_sect):
      assert("Different number of bonds")
   
    #"Hashing atom numbers to atom names"
    for j, (old_atoms, new_atoms) in enumerate(zip(BND_old_sect[0:section_elements], BND_new_sect[0:section_elements])):
      for i, (old_atom, new_atom) in enumerate(zip(old_atoms, new_atoms)):
        old_atom = int(old_atom)
        new_atom = int(new_atom)
        # new_atom = int(new_atom)
        # if(old_atom<0):
        #   old_atom = abs(old_atom)
        # if(new_atom<0):
        #   new_atom = abs(new_atom)
        BND_old_sect[j][i] = self.old_names[(old_atom)-1].strip('_')
        BND_new_sect[j][i] = self.new_names[(new_atom)-1].strip('_')


    BND_old_sect = zip(*BND_old_sect)
    BND_new_sect = zip(*BND_new_sect)

    return BND_old_sect, BND_new_sect


def _file_len(fname):
  with open(fname) as f:
    for i, l in enumerate(f):
      pass
  return i + 1

