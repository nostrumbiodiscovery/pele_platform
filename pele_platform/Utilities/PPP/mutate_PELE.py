import argparse
import os
import sys
from multiprocessing import Pool, cpu_count
import subprocess
import PPP.mut_prep4pele as mut
import PPP.enviroment_parameters as env
import Utilities.template_builder as tb
import Utilities.helpers as hp
import Utilities.best_structs as bs
import shutil
import AdaptivePELE.adaptiveSampling as ad

PELE_BIN = os.path.join(env.pele_folder_path, "bin/Pele_mpi")
DIR = os.path.dirname(os.path.realpath(__file__))
ADAPTIVE_CONTROL_FILE = os.path.join(DIR, "Template/adaptive.conf")
PELE_CONTROL_FILE = os.path.join(DIR, "Template/pele.conf")
CMD = "python -m AdaptivePELE.adaptiveSampling {}"


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_system", help="Input to apply the mutations to")
    parser.add_argument("file", help="File to read mutations from. Each line needs to follow the next format:  mut1.pdb ARG 164 ALA", type=str)
    parser.add_argument("--box", nargs="+", help="Center of the box")
    parser.add_argument("--cpus", type=int, help="Cpus to use", default=1)
    args = parser.parse_args()
    return args.input_system, args.file, args.box, args.cpus

def main(system, mutations_file, box, cpus):
    mutations = retrieve_mutation(mutations_file)
    for mutants in mutations:
        mutants_for_pele = []
        for mutation_info in mutants:
            print(mutation_info.split(" ",1))
            output, mutation = mutation_info.split(" ",1)
            ini_res, resnum, chain, end_res = mutation.split()
            mutation = {'fin_resname': end_res, 'resnum': resnum, 'chain': "", 'ini_resname': ini_res}
            mutants_for_pele.append(mut.main(system, 2.5, mutation=[mutation], output_pdb=[output])[0])
        folders = [ os.path.abspath(os.path.splitext(mutant)[0]) for mutant in mutants_for_pele ]
        pool = Pool(cpus)
        pool.map(launch_mutation_PELE, mutants_for_pele)
        pool.close()
        system = bs.main(path=folders, criteria="Binding Energy", n_structs=1, sort_order="min", out_freq=1, output="Mutation")




def retrieve_mutation(mutations_file):
    with open(mutations_file, "r") as f:
        line = f.readline()
        mutations = []
        mutation = []
        while line:
            if "MUTATION" not in line:
                if line.strip():
                    mutation.append(line.strip("\n").strip().strip("\t"))
            else:
                if mutation:
                    mutations.append(mutation)
                mutation = []
            line = f.readline()
        mutations.append(mutation)
    return mutations



def launch_mutation_PELE(mutant):
    root_name = os.path.abspath(os.path.splitext(mutant)[0])
    ad_conf = prepare_control_file(root_name, mutant, box, adaptive=True)
    _ = prepare_control_file(root_name, mutant, box, adaptive=False)
    subprocess.call(CMD.format(ad_conf).split())

def prepare_control_file(root_name, mutant, box, adaptive=False):
    conf_name = root_name + "_adaptive.conf" if adaptive else root_name + "_pele.conf"
    KEYWORDS = {"OUTPUT": root_name, "STRUCT": os.path.abspath(mutant), "PELE": root_name + "_pele.conf"} if adaptive else {"BOX" : ",".join(box)}
    template_file = ADAPTIVE_CONTROL_FILE if adaptive else  PELE_CONTROL_FILE
    shutil.copy(template_file, conf_name)
    tb.TemplateBuilder(conf_name, KEYWORDS)
    return conf_name


if __name__ == "__main__":
    system, mutations, box, cpus = parse_args()
    main(system, mutations, box, cpus)
