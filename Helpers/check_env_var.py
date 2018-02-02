import os
import sys
import constants

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


def check_dependencies():

        os.environ["SCHRODINGER"] = constants.SCHRODINGER
        os.environ["PELE"] = constants.PELE
        os.environ["PATH"] = "{}:{}".format(os.environ["PATH"], constants.MPIRUN)

        sys.path.insert(0, os.path.join(os.environ["SCHRODINGER"], "internal/lib/python2.7/site-packages/"))       
        sys.path.insert(0, constants.ADAPTIVE)        
        """
        try:
            os.environ["LD_LIBRARY_PATH"] = "{1}:{0}".format(os.environ["LD_LIBRARY_PATH"], os.path.join(os.environ["SCHRODINGER"],"mmshare-v4.0/lib/Linux-x86_64/"))
        except KeyError:
            os.environ["LD_LIBRARY_PATH"] = os.path.join(os.environ["SCHRODINGER"],"mmshare-v4.0/lib/Linux-x86_64/")
        """
        try:
		os.environ["SCHRODINGER"]
	except KeyError:
		raise("SCHRODINGER IS NOT EXPORTED. export SCHRODINGER='/path/to/schrodinger/")

	try:
		os.environ["PELE"]
		os.environ["PATH"]="{}:{}".format(
			os.environ["PATH"],
			os.path.join(os.environ["PELE"], "bin"))
	except KeyError:
		raise("PELE IS NOT EXPORTED. Export PELE='/path/to/pele/")

	if not find_executable("mpirun"):
		raise ValueError("Could not find executable in path. Set export PATH=/path/to/mpirun_binary/:$PATH")

	if not find_executable("Pele_mpi") and not find_executable("PELE-1.5_mpi"):
		raise ValueError("Could not find executable in path. set export PATH=/path/to/pele_binary/:$PATH")
	


