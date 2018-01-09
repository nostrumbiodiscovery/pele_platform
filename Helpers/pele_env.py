import sys
import os
import shutil
import warnings



class Pele_env_Builder(object):
		"""
			Base class wher the needed pele environment
			is build by creating folders and files
		"""

		def __init__(self, input, forcefield, template, rotamers, pele_dir):
			self.input = input
			self.forcefield = forcefield
			self.template = template
			self.rotamers = rotamers
			self.pele_dir = pele_dir
			self.templates = os.path.abspath(os.path.join(os.path.dirname(os.path.dirname(__file__)), "PeleTemplates"))

		def folder_levels(self):
			"""
				Create pele folders
			"""

			folders_to_create = ["",
								 "DataLocal/Templates/OPLS2005/HeteroAtoms/",
								 "DataLocal/Templates/AMBER99sb/HeteroAtoms/",
								 "DataLocal/Templates/AMBER99sbBSC0/HeteroAtoms/",
								 "DataLocal/LigandRotamerLibs",
								 "results",
								 "output",
								 "output_adaptive_short",
								 "output_clustering",
								 "output_adaptive_long",
								 ]
			
			for folder in folders_to_create:
				self.create_dir(self.pele_dir, folder)
			

		def file_dist(self):
			"""
				Copy control rotamer
				and template files
			"""

			#paths
			self.template_dir = os.path.join(
				self.pele_dir, "DataLocal/Templates/{}/HeteroAtoms/".format(self.forcefield))
			self.rotamers_dir = os.path.join(
				self.pele_dir, "DataLocal/LigandRotamerLibs")

			#actions
			shutil.copy(self.input, os.path.join(self.pele_dir, "complex.pdb"))
			with cd(os.path.abspath(self.templates)):
					self.copy("adaptive_short.conf", self.pele_dir)
					self.copy("pele.conf", self.pele_dir)
					self.copy("adaptive_long.conf", self.pele_dir)
			
			shutil.move(self.template, self.template_dir)
			shutil.move(self.rotamers, self.rotamers_dir)

			adaptive_short_temp = os.path.join(self.pele_dir,"adaptive_short.conf")
			pele_temp = os.path.join(self.pele_dir,"pele.conf")
			adaptive_long_temp = os.path.join(self.pele_dir,"adaptive_long.conf")

			return adaptive_short_temp, pele_temp, adaptive_long_temp




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
 					warnings.warn("Directory {} already exists.".format(path),RuntimeWarning)
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

class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)


def set_pele_env(system,  forcefield, template, rotamers_file, pele_dir):
	pele_env = Pele_env_Builder(system,  forcefield, template, rotamers_file, pele_dir)
	pele_env.folder_levels()
	adaptive_short_temp, pele_temp, adaptive_long_temp = pele_env.file_dist()
	return adaptive_short_temp, pele_temp, adaptive_long_temp
