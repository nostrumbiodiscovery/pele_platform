import os, errno

def silentremove(*args, **kwargs):
    for files in args:
        for filename in files:
    	    try:
        	os.remove(filename)
    	    except OSError as e: # this would be "except OSError, e:" before Python 2.6
        	if e.errno != errno.ENOENT: # errno.ENOENT = no such file or directory
            	    raise # re-raise exception if a different error occurred


class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)

