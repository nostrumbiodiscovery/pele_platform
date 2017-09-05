import pytest
import sys
import os
import subprocess
sys.path.insert(1, os.path.join(sys.path[0], '..'))
#from ciphers.DESchiper import convert_to_binary, DES, XOR

SCHRODINGER_PYTHON = "/opt/schrodinger2016-4/utilities/python"

@pytest.mark.parametrize("argument, expected", [
                         ('101101', '011011'),
                         ('111101', '111011'),
                         ('101000', '100010'),
                         ])
def test_keyrotation(argument, expected):
    subprocess.call([SCHRODINGER_PYTHON, "main.py", argument])
    key_rotated = DES_cipher.key_rotation(key)
    assert key_rotated == bt(expected)