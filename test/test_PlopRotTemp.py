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
BONDS = [[0, 1], [1, 2], [1, 3], [3, 4], [3, 8], [4, 5], [4, 13], [5, 6], 
        [5, 14], [6, 7], [6, 15], [7, 8], [7, 16], [8, 9], [9, 10], [10, 11],
        [10, 12], [12, 17], [12, 18], [12, 19]]
PARENT_RESULT = [0, 1, 2, 2, 4, 5, 6, 7, 4, 9, 10, 11, 11, 5, 6, 7, 8, 13, 13, 13]



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
# def test_build_template(MAE_FILE, root, OPLS, hetgrp_opt, old_name, new_name, expected):
#     [template_file, output_template_file, mae_file_hetgrp_ffgen, files, resname] = pl.build_template(MAE_FILE, root, OPLS, hetgrp_opt, old_name, new_name)
#     assert template_file == expected




@pytest.mark.parametrize("mae_file, template_file, mae_expected, template_expected", [
                         (os.path.join(TEST_PATH, 'ain.mae'), 'ain.hetgrp_ffgen' , MAE_CONVERSION, TEMPLATE_CONVERSION),
                         ])
def test_MatchTempMaeAtoms(mae_file, template_file, mae_expected, template_expected):
    [mae2temp, temp2mae] = pl.MatchTempMaeAtoms(mae_file, template_file)
    assert mae2temp == template_expected
    assert temp2mae == mae_expected

@pytest.mark.parametrize("mae_file, expected", [
                         (os.path.join(TEST_PATH, 'ain.mae'), ''),
                         (os.path.join(TEST_PATH, 'ain_repited.mae'), Exception),
                         ])
def test_check_repite_names(mae_file, expected):
    atomnames = pl.find_names_in_mae(mae_file)
    if(expected == Exception):
        with pytest.raises(expected):
            pl.check_repite_names(atomnames)
    else:
        pass



@pytest.mark.parametrize("rotamer_library", [
                         (os.path.join(TEST_PATH,'ain_vdw')),
                         ])
def test_check_replace_vdwr(rotamer_library):
    pl.replace_vdwr_from_library(rotamer_library)
    radius_vdw_info, start_index, end_index = pl.parse_nonbonded(rotamer_library)
    for i, rdw_line in enumerate(radius_vdw_info):
        NBOND_info = rdw_line.split()
        rdw = float(NBOND_info[1])/2.0
        if(rdw== 0):
            assert 0
 
@pytest.mark.parametrize("bonds, n_atoms, expected", [
                         (BONDS, 20, PARENT_RESULT),
                         ])
def test_search_parent_atom(bonds, n_atoms, expected):
    parents = pl.search_parent_atom(bonds, n_atoms)
    assert parents == expected
    



