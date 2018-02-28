import os
import shutil
import warnings


def set_pele_env(folders, files, forcefield, pele_dir):
    pele_env = Pele_env_Builder(folders, files, forcefield, pele_dir)
    pele_env.folder_levels()
    pele_env.file_dist()


class Pele_env_Builder(object):
    """
        Base class wher the needed pele environment
        is build by creating folders and files
    """

    def __init__(self,folders, files, forcefield, pele_dir):
        self.folders = folders
        self.files = files
        self.forcefield = forcefield
        self.pele_dir = pele_dir
        self.templates = os.path.abspath(os.path.join(os.path.dirname(os.path.dirname(__file__)), "PeleTemplates"))

    def folder_levels(self):
        """
            Create pele folders
        """

        for folder in self.folders:
            self.create_dir(self.pele_dir, folder)

    def file_dist(self):
        """
            Copy templates
        """

        # Actions
        for file in self.files:
            self.copy(file, self.pele_dir)


    def create_dir(self, base_dir, extension=None):
        """
            Class Method to manage
            directory creation only if that
            ones doesn't exist

            Location:
                base_dir+extension
                or base_dir if extension is None
        """
        if extension:
            path = os.path.join(base_dir, extension)
            if os.path.isdir(path):
                warnings.warn("Directory {} already exists.".format(path), RuntimeWarning)
            else:
                os.makedirs(path)
        else:
            if os.path.isdir(base_dir):
                warnings.warn("Directory {} already exists.".format(base_dir), RuntimeWarning)
            else:
                os.makedirs(base_dir)

    def copy(self, standard, destination, user=None):
        if user:
            shutil.copy(user, os.path.join(self.pele_dir, standard))
        else:
            shutil.copy(standard, self.pele_dir)
        return os.path.join(self.pele_dir, standard)


def is_repited(pele_dir):

    original_dir = None
    split_dir = pele_dir.split("_")
    for chunk in split_dir:
        if chunk != "Pele":
            if original_dir:
                original_dir = "{}_{}".format(original_dir, chunk)
            else:
                original_dir = chunk
        else:
            break
    if split_dir[-1].isdigit():
        i = split_dir[-1]
        i = int(i) + 1 
    else:
        i = 1
    if os.path.isdir(pele_dir):
		new_pele_dir = "{}_Pele_{}".format(original_dir, i)
		new_pele_dir = is_repited(new_pele_dir)
		return new_pele_dir
    else:
		return pele_dir

def is_last(pele_dir):

    original_dir = None
    split_dir = pele_dir.split("_")
    for chunk in split_dir:
		if chunk != "Pele":
			if original_dir:
 				original_dir = "{}_{}".format(original_dir, chunk)
			else:
				original_dir = chunk
		else:
			break
    if split_dir[-1].isdigit():
        i = split_dir[-1]
        i = int(i) + 1 
    else:
		i = 1 

    if os.path.isdir(pele_dir):
            new_pele_dir = "{}_Pele_{}".format(original_dir, i)
            if not os.path.isdir(new_pele_dir):
                return pele_dir
            else:
			    new_pele_dir = is_last(new_pele_dir)
			    return new_pele_dir
    else:
        return pele_dir

