import pytest
import sys
import os
import subprocess
sys.path.insert(1, os.path.join(sys.path[0], '..'))
import PlpRotTemp.PlopRotTemp as pl
#from ciphers.DESchiper import convert_to_binary, DES, XOR

TEST_PATH = os.path.join(os.path.dirname(__file__), 'data')
MAE_FILE = 'ain.mae'
ROOT = 'ain'
MAE_CONVERSION = [0,1,2,3,4,13,5,14,6,15,7,16,8,9,10,11,12,17,18,19]
TEMPLATE_CONVERSION = [0,1,2,3,4,6,8,10,12,13,14,15,16,5,7,9,11,17,18,19]

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
    [template_file, output_template_file, mae_file_hetgrp_ffgen, files, resname] = pl.build_template(MAE_FILE, root, OPLS, hetgrp_opt, old_name, new_name)
    assert template_file == expected

"""
@pytest.mark.parametrize("mae_min_file, user_fixed_bonds, log_file, use_rings, \
                          use_mult_lib, user_core_atom, user_tors, back_tors, \
                          max_tors, R_group_root_atom_name", 
                          [
                         (os.path.join(TEST_PATH, 'ain.mae'), 'ain', '2005', '', '', '', 'ain.hetgrp_ffgen'),
                         (os.path.join(TEST_PATH,'MI4.mae'), 'MI4', '2005', '', '', '', 'mi4.hetgrp_ffgen'),
                         ])

def test_FindCore(mae_min_file, user_fixed_bonds, log_file, use_rings, \
                 use_mult_lib, user_core_atom, user_tors, back_tors, max_tors, R_group_root_atom_name):
    [mae_num, parent, rank, tors, use_rings, group, back_tors, tors_ring_num] = \
        pl.FindCore(mae_min_file, user_fixed_bonds, log_file, use_rings, \
                 use_mult_lib, user_core_atom, user_tors, back_tors, max_tors, R_group_root_atom_name)
    assert template_file == expected
"""


@pytest.mark.parametrize("mae_file, template_file, mae_expected, template_expected", [
                         (os.path.join(TEST_PATH, 'ain.mae'), 'ain.hetgrp_ffgen' , MAE_CONVERSION, TEMPLATE_CONVERSION),
                         ])
def test_MatchTempMaeAtoms(mae_file, template_file, mae_expected, template_expected):
    [mae2temp, temp2mae] = pl.MatchTempMaeAtoms(mae_file, template_file)
    assert mae2temp == template_expected
    assert temp2mae == mae_expected

"""###############testttt
@pytest.mark.parametrize("mae_file, template_file, mae_expected, template_expected", [
                         (os.path.join(TEST_PATH, 'ain.mae'), 'ain.hetgrp_ffgen' , MAE_CONVERSION, TEMPLATE_CONVERSION),
                         ])
def test_make_libraries(mae_file, template_file, mae_expected, template_expected):
    rotamers_file = pl.make_libraries(resname, conf_file, root, names, zmat_atoms,
                                      group, use_rings, use_mult_lib,
                                      output_template_file, gridres, debug)
    assert mae2temp == template_expected
    assert temp2mae == mae_expected
"""


