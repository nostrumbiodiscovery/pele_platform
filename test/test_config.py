import sys
import os
import re
import subprocess

MAE_FILE = 'ain.mae'
ROOT = 'ain'
MAE_CONVERSION = [0,1,2,3,4,13,5,14,6,15,7,16,8,9,10,11,12,17,18,19]
TEMPLATE_CONVERSION = [0,1,2,3,4,6,8,10,12,13,14,15,16,5,7,9,11,17,18,19]
BONDS = [[0, 1], [1, 2], [1, 3], [3, 4], [3, 8], [4, 5], [4, 13], [5, 6], 
        [5, 14], [6, 7], [6, 15], [7, 8], [7, 16], [8, 9], [9, 10], [10, 11],
        [10, 12], [12, 17], [12, 18], [12, 19]]
PARENT_RESULT = [0, 1, 2, 2, 4, 5, 6, 7, 4, 9, 10, 11, 11, 5, 6, 7, 8, 13, 13, 13]
REPOS_PATH = os.path.abspath(os.path.join(__file__ ,"../.."))
TEST_PATH = os.path.join(os.path.dirname(__file__), 'data')
MAIN_PATH =  os.path.join(REPOS_PATH, 'PlpRotTemp/main.py')
OLD_MAIN_PATH = '/home/dani/repos/presentation/PlpRotTemp/main_16_9_17.py'
try:
    PYTHON_PATH = os.path.join(os.environ['SCHRODINGER'] + "/utilities/python")
except KeyError:
    print("Set SCHRODINGER environment variable path")
TRESHOLD= 0.15



def run_versions(input_file):
  try:
    os.remove("lig")
    os.remove("LIG")
  except OSError:
    pass
  res_name = "LIG"
  NAMES = ['atom1', 'atom2', 'spring', 'eq_dist']
  subprocess.call([PYTHON_PATH, MAIN_PATH, os.path.join(TEST_PATH, input_file)])
  subprocess.call([PYTHON_PATH, OLD_MAIN_PATH, os.path.join(TEST_PATH, input_file)])
  new_sections = parse_template(os.path.join(REPOS_PATH,res_name.upper()), res_name.upper())
  old_sections = parse_template(os.path.join(REPOS_PATH,res_name.lower()), res_name.upper())
  return new_sections, old_sections

def parse_template(file_path, input_file):
  
  res_name = input_file.split('.')[0]

  starts = [False, False, False, False, False, False]
  ends = [False, False, False, False, False, False]
  KEYWORDS = [res_name, 'NBON', 'BOND', 'PHI', 'IPHI']
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
  

  with open(file_path, 'r') as f:
    
    line = f.readline() #Skip Commentaris
    line = f.readline() #Skip Commentaris
    line = f.readline() #Skip Commentaris
    for i, (keyword, item, index)  in enumerate(zip(
        KEYWORDS, ITEMS, INDEXES)):
        line = re.sub(' +',' ',line)
        line = line.strip('\n').strip()
        while not ends[i]:
            # print(line)
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
    # print(ITEMS)
    return ITEMS


