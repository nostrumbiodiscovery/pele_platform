import pytest
import sys
import os
import subprocess
sys.path.insert(1, os.path.join(sys.path[0], '..'))
import PlpRotTemp.PlopRotTemp as pl
#from ciphers.DESchiper import convert_to_binary, DES, XOR

SCHRODINGER_PYTHON = "/opt/schrodinger2016-4/utilities/python"
MAE_FILE = 'ain.mae'
ROOT = 'ain'
TEST_PATH = os.path.dirname(__file__)

@pytest.mark.parametrize("argument, expected", [
                         ('ain.mae', 'ain'),
                         ('~/ain.mae', 'ain'),
                         ('/opt/schrodinger/ain.mae', 'ain'),
                         ])
def test_get_root(argument, expected):
    root = pl.get_root_path(argument)
    assert root == expected

@pytest.mark.parametrize("MAE_FILE, root, OPLS, hetgrp_opt, old_name, new_name, expected", [
                         (os.path.join(TEST_PATH, 'ain.mae'), 'ain', '2005', '', '', '', 'ain.hetgrp_ffgen'),
                         (os.path.join(TEST_PATH,'MI4.mae'), 'MI4', '2005', '', '', '', 'mi4.hetgrp_ffgen'),
                         ])
def test_build_template(MAE_FILE, root, OPLS, hetgrp_opt, old_name, new_name, expected):
    print(MAE_FILE)
    [template_file, output_template_file, mae_file_hetgrp_ffgen, files, resname] = pl.build_template(MAE_FILE, root, OPLS, hetgrp_opt, old_name, new_name)
    assert template_file == expected

