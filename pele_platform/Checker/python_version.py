import sys
import pele_platform.Errors.custom_errors as ce

def check_python_version():
    python_version = sys.version_info
    if sys.version_info < (3, 0):
        raise ce.OldPythonVersion("Pele Platform requieres Python>3.x. Used is {}.{}".format(
        python_version.major, python_version.minor))
