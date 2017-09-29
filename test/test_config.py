import sys
import os
import re



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


