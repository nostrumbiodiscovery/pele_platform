import sys
from tmp_helpers import Helper



class ChargeHandler(Helper):

	def __init__(self, file, num_atoms):
		self.file = file
		self.num_atoms = num_atoms
		

	def get_charges(self):


		lines = self.preproces_file_lines(self.file)

		for line in lines:
			if line: print(line)
		charges = [line.replace(',','.') for line in lines if line]

		if len(charges)!= self.num_atoms: 
			raise ValueError("Not the same number of Charges and pdb Atoms")
		else:
			return charges




	
