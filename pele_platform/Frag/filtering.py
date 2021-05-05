import argparse
import glob
import os
import yaml
import time
import pickle
import sys

from drug_learning.two_dimensions.Input import fingerprints as fp
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import QED
from rdkit.Chem import RDConfig
from yaml import load, Loader
from collections import ChainMap
from multiprocessing import Pool
from functools import partial


class Library:

    def __init__(self, path, init_mol, filters=None, save=False):
        self.path = path  # not using absolute path since we defined path on parse_args()
        self.filters = filters
        self.errors = 0
        self.mol_num = 0
        self.filtered_final = []
        # self.sd_files=self._retrieve_files()
        if init_mol:
            mol_descriptors = {}
            frag = Fragment(init_mol)
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
        mor = fp.MorganFP()
        structures = mor.fit(input_sdf)
        features = mor.transform()
        return features

    def main(self, files_sdf):
        fingerprints = {}
        for file_sdf in files_sdf:
            try:
                counter = 0
                fingerprint = self.generate_fingerprint(file_sdf)
                mols = Chem.SDMolSupplier(file_sdf, removeHs=False)
                name = os.path.basename(file_sdf).rsplit(".")[0]
            except:
                continue
            for i in range(len(mols)):
                fingerprints[mols[i]] = fingerprint[i]
                self.molecule = mols[i]
                self.fragments_dum = [Fragment(mols[i])]
                self.filters_f()
                self.save()
            print("\ (•◡•) / File %s finished \ (•◡•) / " % name)
            print(".。・゜・。..。・゜・。..。・゜・。..。・゜・。..。・゜・。..。・゜・。..。・゜・。.")

    def filters_f(self):
        self.filtered_fragments = self._apply_filters()

    def save(self):
        for mol in self.filtered_fragments:
            self.filtered_final.append(mol)

    def _retrieve_files(self):  # deleted argument path since we can access it within method
        sdf_path = os.path.join(self.path, "*.sdf")
        sd_files = glob.glob(sdf_path)

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
                # import pdb
                # pdb.set_trace()
                filter_pass = []
                # print('mw:',float(self.parsed_filters['mw'][0])<frag.mw<float(self.parsed_filters['mw'][1]))
                # print('logP:,',abs(float(self.parsed_filters['logP'][0]))<abs(frag.logP)<abs(float(self.parsed_filters['logP'][1])))
                # print('hbd:',int(self.parsed_filters['hbd'][0])<frag.hbd<int(self.parsed_filters['hbd'][1]))
                # print('hba:',int(self.parsed_filters['hba'][0])<frag.hba<int(self.parsed_filters['hba'][1]))
                # print('psa:',int(self.parsed_filters['psa'][0])<frag.psa<int(self.parsed_filters['psa'][1]))
                # print('rotb:',int(self.parsed_filters['rotb'][0])<frag.rotb<int(self.parsed_filters['rotb'][1]))

                if 'mw' in self.parsed_filters:
                    # print(float(self.parsed_filters['mw'][0]),frag.mw,float(self.parsed_filters['mw'][1]))
                    # sys.exit()
                    if float(self.parsed_filters['mw'][0]) < frag.mw < float(self.parsed_filters['mw'][1]):
                        filter_pass.append(True)
                    else:
                        filter_pass.append(False)
                        #print('mw:', float(self.parsed_filters['mw'][0]) , frag.mw , float(self.parsed_filters['mw'][1]))
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
                        #print('hbd:', int(self.parsed_filters['hbd'][0]) , frag.hbd , int(self.parsed_filters['hbd'][1]))
                if 'hba' in self.parsed_filters:
                    if int(self.parsed_filters['hba'][0]) < frag.hba < int(self.parsed_filters['hba'][1]):
                        filter_pass.append(True)
                    else:
                        filter_pass.append(False)
                        #print('hba:', int(self.parsed_filters['hba'][0]) , frag.hba , int(self.parsed_filters['hba'][1]))
                if 'psa' in self.parsed_filters:
                    if int(self.parsed_filters['psa'][0]) < frag.psa < int(self.parsed_filters['psa'][1]):
                        filter_pass.append(True)
                    else:
                        filter_pass.append(False)
                        #print('psa:', int(self.parsed_filters['psa'][0]) , frag.psa , int(self.parsed_filters['psa'][1]))
                if 'rotb' in self.parsed_filters:
                    if int(self.parsed_filters['rotb'][0]) <= frag.rotb <= int(self.parsed_filters['rotb'][1]):
                        filter_pass.append(True)
                    else:
                        filter_pass.append(False)
                        #print('rotb:', int(self.parsed_filters['rotb'][0]) , frag.rotb , int(self.parsed_filters['rotb'][1]))
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


def main(ligand, path, filters):
    print("INITIALIZING FILTERING")
    mol = Chem.MolFromPDBFile(ligand, removeHs=False)
    Chem.MolToPDBFile(mol, "original.pdb")

    lib = Library(path, mol, filters)
    return lib.filtered_final
