# -*- coding: utf-8 -*-

import sys
from tmp_helpers import Helper



class ChargeHandler(Helper):

	"""
		Handler class in charge of retrieving
		atom charges from file to TemplateBuilder
	"""

	def __init__(self, file, num_atoms):
		self.file = file
		self.num_atoms = num_atoms
		

	def get_charges(self):
		"""
			Parse self.file for charges 
			and retrieve them if the 
			nº of them is equal to the
			number of atoms in the ligand.
		"""

		lines = self.preproces_file_lines(self.file)

		charges = [line.replace(',','.') for line in lines if line]

		if len(charges)!= self.num_atoms: 
			raise ValueError("Not the same number of Charges and pdb Atoms")
		else:
			return charges




	
