import os
import glob
import shutil
import pele_platform.constants.constants as cs
import pele_platform.main as main

test_path = os.path.join(cs.DIR, "Examples")


FRAG_ARGS = os.path.join(test_path, "frag/input.yaml")
FRAG_SIM_ARGS = os.path.join(test_path, "frag/input_sim.yaml")
FRAG_CORE_ARGS = os.path.join(test_path, "frag/input_core.yaml")
FRAG_AI_ARGS = os.path.join(test_path, "frag/input_ai.yaml")


def test_frag(ext_args=FRAG_ARGS, output="1w7h_preparation_structure_2w_aminoC1N1"):
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

def test_frag_ai(ext_args=FRAG_AI_ARGS, output="round*/"):
    folders = glob.glob(output)
    for folder in folders:
        if os.path.exists(folder):
            shutil.rmtree(folder)
    job = main.run_platform(ext_args)
