import os
import re
import logging

def silentremove(*args, **kwargs):
    for files in args:
        for filename in files:
            try:
                os.remove(filename)
            except OSError:
                pass


class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)


def preproces_lines(lines):
    for i, line in enumerate(lines):
        line = re.sub(' +', ' ', line)
        line = line.strip('\n').strip().split()
        lines[i] = line
    return lines


def set_logger(pele_dir, residue):
	# Logging definition block
    log_name = os.path.abspath("{}.log".format(residue))
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter("%(asctime)s:%(levelname)s:%(message)s")
    file_handler = logging.FileHandler(log_name, mode='w')
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    return logger, log_name

