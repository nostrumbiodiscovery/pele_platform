import os
import sys
import glob
from MSM_PELE import constants


def find_executable(executable, path=None):
    """
    Find if 'executable' can be run. Looks for it in 'path'
    (string that lists directories separated by 'os.pathsep';
    defaults to os.environ['PATH']). Checks for all executable
    extensions. Returns full path or None if no command is found.
    """
    if path is None:
        path = os.environ['PATH']
    paths = path.split(os.pathsep)
    extlist = ['']
    if os.name == 'os2':
        (base, ext) = os.path.splitext(executable)
        # executable files on OS/2 can have an arbitrary extension, but
        # .exe is automatically appended if no dot is present in the name
        if not ext:
            executable = executable + ".exe"
    elif sys.platform == 'win32':
        pathext = os.environ['PATHEXT'].lower().split(os.pathsep)
        (base, ext) = os.path.splitext(executable)
        if ext.lower() not in pathext:
            extlist = pathext
    for ext in extlist:
        execname = executable + ext
        if os.path.isfile(execname):
            return execname
        else:
            for p in paths:
                f = os.path.join(p, execname)
                if os.path.isfile(f):
                    return f
    else:
        return None

def patch_environ():
    """

    Patch current environment variables so Schrodinger can start up and we can import its modules.
    Be warned that calling this function WILL restart your interpreter. Otherwise, Python
    won't catch the new LD_LIBRARY_PATH (or platform equivalent) and Schrodinger won't find its
    libraries during import.

    """

    #Check whether patch_environ was already run
    if "SCHRODINGER" in os.environ:
        return
    else:
        os.environ["SCHRODINGER"] = constants.SCHRODINGER

    #Find schrodinger libraries
    schrodinger_libs_pattern = os.path.join(constants.SCHRODINGER, "mmshare*/lib/Linux-x86_64/")
    schrodinger_libs = glob.glob(schrodinger_libs_pattern)
    schrodinger_libs.append(os.path.join(constants.SCHRODINGER, "internal/lib/ssl"))

    #Update LD_LIBRARY_PATH
    if 'LD_LIBRARY_PATH' in os.environ:
        os.environ['LD_LIBRARY_PATH'] = "{}:{}".format(":".join(schrodinger_libs), os.environ['LD_LIBRARY_PATH'])
    else:
        os.environ['LD_LIBRARY_PATH'] = ":".join(schrodinger_libs)

    #Relunch shell
    os.execve(sys.executable, [sys.executable] + sys.argv, os.environ)


def check_dependencies():
        #Update libraries
        patch_environ()

        #Update env_variables
        if "PYTHONPATH" in os.environ:
            os.environ["PYTHONPATH"] = "{}:{}".format(constants.DIR, os.environ["PYTHONPATH"])
        else:
             os.environ["PYTHONPATH"] = constants.DIR

        os.environ["SCHRODINGER"] = constants.SCHRODINGER
        os.environ["PELE"] = constants.PELE

        # Provisonal line, may be necessary for old schrodinger versions
        if constants.MMSHARE is not None:
            os.environ["MMSHARE_EXEC"] = constants.MMSHARE

        sys.path.append(os.path.join(os.environ["SCHRODINGER"], "internal/lib/python2.7/site-packages/"))

        try:
            os.environ["PATH"] = "{}:{}".format(os.environ["PATH"], constants.MPIRUN)
        except ValueError:
            os.environ["PATH"] = constants.MPIRUN

        #Check dependencies
        try:
            os.environ["SCHRODINGER"]
        except KeyError:
            raise("Change SCHRODINGER path in constants.py module")

        try:
            pele_bin_path = os.path.join(os.environ["PELE"], "bin")
            os.environ["PATH"] = "{}:{}".format(os.environ["PATH"], pele_bin_path)
        except KeyError:
            raise("Change PELE path in constants.py module")

        if not find_executable("mpirun"):
            raise ValueError("Change mpirun path in constants.py module")

        if not find_executable("Pele_mpi") and not find_executable("PELE-1.5_mpi"):
            raise ValueError("Change Pele path in constants.py module")
