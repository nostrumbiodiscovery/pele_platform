import os
import glob
import shutil
import pele_platform.constants.constants as cs
import pele_platform.main as main
from tests import test_adaptive as td

test_path = os.path.join(cs.DIR, "Examples")


FRAG_ARGS = os.path.join(test_path, "frag/input.yaml")
FRAG_SIM_ARGS = os.path.join(test_path, "frag/input_sim.yaml")
FRAG_CORE_ARGS = os.path.join(test_path, "frag/input_core.yaml")
FLAGS_ARGS = os.path.join(test_path, "frag/input_flags.yaml")
FRAG_JOINER_ARGS = os.path.join(test_path, "frag/sdf_joiner/*.yml")


def test_frag_standard(ext_args=FRAG_ARGS, output="1w7h_preparation_structure_2w_aminoC1N1"):
    if os.path.exists(output):
        shutil.rmtree(output)
    job = main.run_platform(ext_args)

def test_frag_sim(ext_args=FRAG_SIM_ARGS, output="1w7h_preparation_structure_2w_aminoC1N1"):
    if os.path.exists(output):
        shutil.rmtree(output)
    job = main.run_platform(ext_args)


def test_frag_core(ext_args=FRAG_CORE_ARGS, output="1w7h_preparation_structure_2w_aminoC1N1"):
    if os.path.exists(output):
        shutil.rmtree(output)
    job = main.run_platform(ext_args)

def test_flags(ext_args=FLAGS_ARGS, output="water_processed_processed_aminoCA1N1"):
    FRAG_FLAGS = ['"seed" : 3000',]
    errors = []
    if os.path.exists(output): shutil.rmtree(output, ignore_errors=True)
    job = main.run_platform(ext_args)
    folder = output
    #if not os.path.exists(os.path.join(folder, "DataLocal/LigandRotamerLibs/STR.rot.assign")) or not os.path.exists(os.path.join(folder, "DataLocal/LigandRotamerLibs/MG.rot.assign")):
        #errors.append("External rotamer flag not working")
    #if not os.path.exists(os.path.join(folder, "DataLocal/Templates/OPLS2005/HeteroAtoms/strz")) or not os.path.exists(os.path.join(folder, "DataLocal/Templates/OPLS2005/HeteroAtoms/mgz")):
        #errors.append("External templates flag not working")
    errors = td.check_file(folder, "control_folder/0_pele_template.conf", td.PELE_VALUES + FRAG_FLAGS, errors)
    errors = td.check_file(folder, "DataLocal/LigandRotamerLibs/SB4.rot.assign", "60", errors)
    assert not errors

def test_sdf_joiner(ext_args=FRAG_JOINER_ARGS):
    files = glob.glob(ext_args)
    for file in files:
        try:
            job = main.run_platform(file)
        except Exception:
            assert False
