import sys
import pele_platform.Errors.custom_errors as ce

def check_python_version():
    #Check python version is > 3.7
    python_version = sys.version_info
    if sys.version_info < (3, 7):
        raise ce.OldPythonVersion("Pele Platform requieres Python>3.6. Used is {}.{}".format(
        python_version.major, python_version.minor))
