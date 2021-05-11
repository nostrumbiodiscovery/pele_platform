import argparse
import os
import yaml
import time
import pickle
import sys
import csv

from drug_learning.two_dimensions.Input import fingerprints as fp
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import QED
from rdkit.Chem import RDConfig
from yaml import load, Loader
from collections import ChainMap
from multiprocessing import Pool
from functools import partial
import openbabel
from openbabel import pybel
from rdkit import DataStructs
from glob import glob


class Library:

    def __init__(self, path, ligand, ligand_path, filters=None, save=False):
        self.path = path  # not using absolute path since we defined path on parse_args()
        self.ligand = ligand
        self.filters = filters
        self.errors = 0
        self.mol_num = 0
        self.filtered_final = []
        #self.init_mol = Chem.MolFromPDBFile(ligand, removeHs=False)
        self.init_mol = ligand
        self.init_mol_path = ligand_path
        self.database_files = glob(path + '/*.csv')[:20]
        # self.sd_files=self._retrieve_files()
        if self.init_mol:
            mol = next(pybel.readfile("pdb", self.init_mol_path))
            finalSDF = pybel.Outputfile("sdf", "input_ligand.sdf", overwrite=True)
            finalSDF.write(mol)
            self.init_mol.fingerprint = self.generate_fingerprint("input_ligand.sdf")
            mol_descriptors = {}
            frag = Fragment(self.init_mol)
            mol_descriptors['mw'] = [frag.mw - 50, frag.mw + 50]
            mol_descriptors['logP'] = [frag.logP - 5, frag.logP + 5]
            mol_descriptors['hbd'] = [frag.hbd - 5, frag.hbd + 5]
            mol_descriptors['hba'] = [frag.hba - 5, frag.hba + 5]
            mol_descriptors['psa'] = [frag.psa - 5, frag.psa + 5]
            mol_descriptors['rotb'] = [frag.rotb - 5, frag.rotb + 5]
            mol_descriptors['arom'] = [frag.arom]
            self.parsed_filters = mol_descriptors
        else:
            self._load_filters()

        self.sd_files = self._retrieve_files()
        self.main(self.sd_files)

    def generate_fingerprint(self, input_sdf):
        try:
            mor = fp.MorganFP()
            structures = mor.fit(input_sdf)
            features = mor.transform()
        except:
            return 'None'
        return features

    def read_fingerprints(self):
        dict = {}
        for f in self.database_files:
            with open(f) as csvfile:
                reader = csv.reader(csvfile, delimiter='\n')
                for row in reader:
                    dict[row[0].split(',')[0]] = row[0].split(',')[1:]
        return dict

    def tanimoto_coefficient(self, input_ligand, database_molecule):
        input_ligand = input_ligand[0].tolist()
        N = len(input_ligand)
        assert N == len(database_molecule)
        v1v2, v1v1, v2v2 = 0., 0., 0.
        for i in range(N):
            v1v2 += input_ligand[i] * int(database_molecule[i])
            v1v1 += input_ligand[i] * input_ligand[i]
            v2v2 += int(database_molecule[i]) * int(database_molecule[i])
        return v1v2 / (v1v1 + v2v2 - v1v2)


    def main(self, files_sdf):
        self.fingerprints = self.read_fingerprints()
        self.names = {}
        self.coefficients = {}
        for file_sdf in files_sdf:
            mols = Chem.SDMolSupplier(file_sdf, removeHs=False)
            name = os.path.basename(file_sdf).rsplit(".")[0]
            for mol in mols:
                if mol:
                    self.names[mol.GetProp("_Name")] = mol
                    self.coefficients[mol.GetProp("_Name")]= self.tanimoto_coefficient(self.init_mol.fingerprint, self.fingerprints[mol.GetProp("_Name")])
                    print("Tanimoto coefficient: %s" % self.coefficients[mol.GetProp("_Name")])

                #         if mol:
        #             Chem.MolToPDBFile(mol, "molecule_database.pdb")
        #             mol_to_sdf = next(pybel.readfile("pdb", "molecule_database.pdb"))
        #             finalSDF = pybel.Outputfile("sdf", "database_molecule.sdf", overwrite=True)
        #             finalSDF.write(mol_to_sdf)
        #             mol.fingerprint = self.generate_fingerprint("database_molecule.sdf")
        #             if mol.fingerprint == 'None':
        #                 continue
        #             mol.tanimoto = self.tanimoto_coefficient(mol.fingerprint, self.init_mol.fingerprint)
        #             coefficients[mol] = mol.tanimoto
        #             print( "Tanimoto coefficient: %s" % mol.tanimoto)
        #             #self.molecule = mol
        #             #self.fragments_dum = [Fragment(mol)]
        #             #self.filters_f()
        #             #self.save()
            print("\ (•◡•) / File %s finished \ (•◡•) / " % name)
            print(".。・゜・。..。・゜・。..。・゜・。..。・゜・。..。・゜・。..。・゜・。..。・゜・。.")
        #top_similar = max(self.coefficients, key=self.coefficients.get)
        top_similar = sorted(self.coefficients, key=self.coefficients.get, reverse=True)[:(int(len(self.coefficients)/4))]
        print('top_similar:', top_similar)
        for molecule in top_similar:
            substructure_results = self.names[molecule].GetSubstructMatches(self.init_mol)
            print('substructure_results:', substructure_results)
        print(len(self.names))

    def filters_f(self):
        self.filtered_fragments = self._apply_filters()

    def save(self):
        for mol in self.filtered_fragments:
            self.filtered_final.append(mol)

    def _retrieve_files(self):  # deleted argument path since we can access it within method
        sdf_path = os.path.join(self.path, "*.sdf")
        sd_files = glob(sdf_path)

        return sd_files

    def _load_filters(self):
        with open(self.filters, 'r') as filters_file:
            yaml = load(filters_file, Loader=Loader)
            filters = yaml["filters"]
        filters_dict = {}
        for k in filters:
            for key, value in filters[k].items():
                try:
                    filters_dict[k].append(value)
                except:
                    filters_dict[k] = []
                    filters_dict[k].append(value)
        return filters_dict

    def _apply_filters(self):
        filtered = []

        for frag in self.fragments_dum:
            if frag.molecule_name:
                filter_pass = []
                if 'mw' in self.parsed_filters:
                    if float(self.parsed_filters['mw'][0]) < frag.mw < float(self.parsed_filters['mw'][1]):
                        filter_pass.append(True)
                    else:
                        filter_pass.append(False)
                if 'logP' in self.parsed_filters:
                    if float(self.parsed_filters['logP'][0]) < frag.logP < float(self.parsed_filters['logP'][1]):
                        filter_pass.append(True)
                    else:
                        filter_pass.append(False)
                if 'hbd' in self.parsed_filters:
                    if int(self.parsed_filters['hbd'][0]) < frag.hbd < int(self.parsed_filters['hbd'][1]):
                        filter_pass.append(True)
                    else:
                        filter_pass.append(False)
                if 'hba' in self.parsed_filters:
                    if int(self.parsed_filters['hba'][0]) < frag.hba < int(self.parsed_filters['hba'][1]):
                        filter_pass.append(True)
                    else:
                        filter_pass.append(False)
                if 'psa' in self.parsed_filters:
                    if int(self.parsed_filters['psa'][0]) < frag.psa < int(self.parsed_filters['psa'][1]):
                        filter_pass.append(True)
                    else:
                        filter_pass.append(False)
                if 'rotb' in self.parsed_filters:
                    if int(self.parsed_filters['rotb'][0]) <= frag.rotb <= int(self.parsed_filters['rotb'][1]):
                        filter_pass.append(True)
                    else:
                        filter_pass.append(False)
                if 'arom' in self.parsed_filters:
                    if frag.arom > 1:
                        filter_pass.append(True)
                    else:
                        filter_pass.append(False)
                if all(filter_pass):
                    filtered.append(frag.molecule)
        return filtered


class Fragment():

    def __init__(self, mol=None, mw=None, logP=None, hba=None, hbd=None, psa=None, rotb=None, arom=None):
        if mol:
            self.molecule = mol if mol else None
            self.molecule_name = mol.GetProp('_Name') if mol.HasProp('_Name') else None
            self._qed = QED.properties(self.molecule)
            self.mw = self._qed.MW
            self.logP = self._qed.ALOGP
            self.hba = self._qed.HBA
            self.hbd = self._qed.HBD
            self.psa = self._qed.PSA
            self.rotb = self._qed.ROTB
            self.arom = self._qed.AROM
        else:
            self.molecule = None
            self.molecule_name = None
            self._qed = None
            self.mw = mw
            self.logP = logP
            self.hba = hba
            self.hbd = hbd
            self.psa = psa
            self.rotb = rotb
            self.arom = arom

    def __str__(self):  # changed to str to ensure readability
        return "Molecule {name}\nMW = {mw}\nlogP = {logp}\n".format(name=self.molecule_name, mw=self.mw, logp=self.logP)


def main(ligand, ligand_path, path, filters):
    print("INITIALIZING FILTERING")
    #mol = Chem.MolFromPDBFile(ligand, removeHs=False)
    #Chem.MolToPDBFile(ligand, "original.pdb")

    lib = Library(path, ligand, ligand_path, filters)
    return lib.filtered_final
