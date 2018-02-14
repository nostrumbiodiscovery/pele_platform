import sys
import os
import shutil
import warnings


class Pele_env_Builder(object):
    """
        Base class wher the needed pele environment
        is build by creating folders and files
    """

    def __init__(self, Input, folders, files, forcefield, template, rotamers, pele_dir):
        self.input = Input
        self.folders = folders
        self.files = files
        self.forcefield = forcefield
        self.template = template
        self.rotamers = rotamers
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
            Copy control rotamer
            and template files
        """
        # Paths
        self.template_dir = os.path.join(
            self.pele_dir, "DataLocal/Templates/{}/HeteroAtoms/".format(self.forcefield))
        self.rotamers_dir = os.path.join(self.pele_dir, "DataLocal/LigandRotamerLibs")

        # Actions
        for file in self.files:
            self.copy(file, self.pele_dir)

        if self.template and self.rotamers:
            shutil.move(self.template, self.template_dir)
            shutil.move(self.rotamers, self.rotamers_dir)

        shutil.move(self.input, self.pele_dir)

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


def set_pele_env(system,  folders, files, forcefield, template, rotamers_file, pele_dir):
    pele_env = Pele_env_Builder(system, folders, files,  forcefield, template, rotamers_file, pele_dir)
    pele_env.folder_levels()
    pele_env.file_dist()

def is_repited(pele_dir):
	if os.path.isdir(pele_dir):
		try:
			original_dir, _, i = pele_dir.rsplit("_", 1)
			i = int(i) + 1 
		except ValueError:
			original_dir, _ = pele_dir.rsplit("_", 1)
			i = 1
		finally:
			new_pele_dir = "{}_Pele_{}".format(original_dir, i)
			new_pele_dir = is_repited(new_pele_dir)
			return new_pele_dir
	else:
		return pele_dir

def is_last(pele_dir):
    if os.path.isdir(pele_dir):
        try:
            original_dir, _, i = pele_dir.r_split("_", 1)
            i = int(i) + 1
        except ValueError:
            original_dir, _ = pele_dir.rsplit("_", 1)
            i = 1
        finally:
            new_pele_dir = "{}_Pele_{}".format(original_dir, i)
            if not os.path.isdir(new_pele_dir):
                return pele_dir
            else:
			    new_pele_dir = is_last(new_pele_dir)
			    return new_pele_dir
    else:
        return pele_dir

