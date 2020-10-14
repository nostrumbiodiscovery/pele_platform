import argparse
from rdkit.Chem import QED
import glob
import os
from rdkit import Chem
import time
from multiprocessing import Pool
import pele_platform.constants.constants as cs

def retrieve_files(path):
    sdf_path = os.path.join(path, "*.sdf")
    sdf_files = glob.glob(sdf_path)
    return sdf_files

def is_donor(props):
    return True if  props.HBD > 0 else False

def is_acceptor(props):
    return True if props.HBA > 0 else False

def is_psa(props, min_psa=0, max_psa=10000):
    return True if min_psa <= props.PSA <= max_psa else False

def is_logp(props, min_logp=-100, max_logp=100):
    return True if min_logp <= props.ALOGP <= max_logp else  False

def is_mw(props, mw=10000):
    return True if props.MW < mw else False
    
def is_rotb(props, rotb=100):
    return True if props.ROTB < rotb else False

def is_aromatic(props):
    return True if props.AROM > 0 else False

def is_toxic(mol):
    
    is_toxic = []

    for index, frag in enumerate(cs.toxic_frags):
        pattern = Chem.MolFromSmarts(frag)

        if frag == "[#16](=[#8])(-[#6])-[#6]": # make sure thioketone is filtered out but not C-SO2-C
            if mol.HasSubstructMatch(pattern) and not mol.HasSubstructMatch(Chem.MolFromSmarts("[#6]-[#16](-[#6])(=[#8])=[#8]")):
                is_toxic.append(mol)
                print("{} has a toxic fragment".format(mol.GetProp("_Name")), "Fragment match:", index, frag)
        else:
            if mol.HasSubstructMatch(pattern):
                is_toxic.append(mol)
                print("{} has a toxic fragment".format(mol.GetProp("_Name")), "Fragment match:", index, frag)

    if is_toxic:
        return True
    else:
        return False

def filter_fragments(f):
    
    restart_folder = "output"
    restart_output = os.path.join(restart_folder, os.path.basename(f))
    print("Running", f)
    if os.path.exists(restart_output):
        return # catch already scanned files
    if not os.path.exists(restart_folder):
        os.mkdir(restart_folder)
    
    writer = Chem.SDWriter(restart_output)
    
    try:
        mols = Chem.SDMolSupplier(f, removeHs=False)
    except OSError:
        return # catch empty files

    for mol in mols:
        if mol:
            atoms = mol.GetAtoms()
            qed = QED.properties(mol)
            if is_mw(qed, 150) and is_donor(qed):
            #if is_donor(qed) and is_acceptor(qed) and is_psa(qed, 0, 10000) and is_logp(qed, -100, 100) and is_mw(qed, 200) and is_rotb(qed, 200) and is_aromatic(qed):
                if not is_toxic(mol):
                    writer.write(mol)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", type=str, required=True, help="Directory with fragment SD files")
    parser.add_argument("--cpus", type=int, help="CPUS to use", default=1)
    parser.add_argument("--save", type=str, required=False, help="SD file name, if you want to save the filtered output.")

    args = parser.parse_args()
    return os.path.abspath(args.dir), args.cpus, args.save

def main():
    sdf_dir, cpus, save = parse_args()
    sdf_list = retrieve_files(sdf_dir)
    with Pool(cpus) as p:
        p.map(filter_fragments, sdf_list)

if __name__ == "__main__":
    main()
