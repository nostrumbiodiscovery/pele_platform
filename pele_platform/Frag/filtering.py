import argparse
import os
import yaml
import time
import pickle
import sys
import csv
import networkx as nx
import shutil
import tempfile
import pandas as pd
import heapq
import numpy as np

from pathlib import Path
from rdkit.Chem import rdMolAlign
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import QED
from rdkit.Chem import RDConfig
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
from rdkit.Chem import rdMolAlign, TemplateAlign
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw
from yaml import load, Loader
from collections import ChainMap
from functools import partial
import openbabel
from openbabel import pybel
from rdkit import DataStructs
from glob import glob
from tqdm import tqdm


class Library:

    def __init__(self, path, ligand, ligand_path, filters=None, save=False):
        self.path = path  # not using absolute path since we defined path on parse_args()
        self.ligand = ligand
        self.filters = filters
        self.errors = 0
        self.mol_num = 0
        self.filtered_final = []
        self.init_mol = ligand
        self.init_mol_path = ligand_path
        self.init_mol_hs = Chem.MolFromPDBFile(self.init_mol_path, removeHs=False)
        self.database_files = glob(path + '/*.csv')
        self.names = {}
        self.coefficients = {}
        self.mols_single_component = {}
        self.substructure_results = {}
        if self.init_mol:
            mol = next(pybel.readfile("pdb", self.init_mol_path))
            finalSDF = pybel.Outputfile("sdf", "input_ligand.sdf", overwrite=True)
            finalSDF.write(mol)
            self.init_mol.fingerprint = Chem.RDKFingerprint(self.init_mol)
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

    def main(self, files_sdf):
        counter = 0
        heap = []
        # Iterate on every sdf file of the database
        for file_sdf in files_sdf:
            mols = Chem.SDMolSupplier(file_sdf, removeHs=True)
            print("COMPUTING TANIMOTO COEFFICIENTS FOR FILE %s: %s" % (file_sdf, counter))
            counter +=1
            for mol in mols:
                if mol:
                    self.molecule = mol
                    self.fragments_dum = [Fragment(mol)]
                    self.filters_f()
                    if self.filters is None:
                        self.filtered_fragments = True
                    if self.filtered_fragments:
                        fp = Chem.RDKFingerprint(mol)
                        coefficient = DataStructs.FingerprintSimilarity(self.init_mol.fingerprint, fp)
                        if coefficient > 0.0:
                            if len(heap) < 100:
                                heapq.heappush(heap, (coefficient, mol.GetProp("_Name")))
                                Chem.MolToPDBFile(mol, "%s.pdb" % (mol.GetProp("_Name")))
                            else:
                                if coefficient > heapq.nsmallest(1, heap)[0][0]:
                                    heapq.heappop(heap)
                                    heapq.heappush(heap, (coefficient, mol.GetProp("_Name")))
                                    try:
                                        Chem.MolToPDBFile(mol, "%s.pdb" % (mol.GetProp("_Name")))
                                    except:
                                        print("FAILED FOR MOL: %s" % (mol.GetProp("_Name")))
                                        continue

        # Only for top similar molecules, perform substructure search
        self.top_molecules = self.get_top_molecule(heap, os.getcwd())
        self.fragment_library = []

        if os.path.exists("./frag_library"):
            shutil.rmtree("./frag_library")
        os.mkdir("./frag_library")
        for ligand in self.top_molecules:
            if os.path.exists(os.path.join("./frag_library", ligand)):
                shutil.rmtree(os.path.join("./frag_library", ligand))
            os.mkdir(os.path.join("./frag_library", ligand))
            with Chem.SDWriter(os.path.join(os.path.join("frag_library/", ligand), ligand + ".sdf")) as w:
                self.ligand_fragment = self.top_molecules[ligand][2]
                self.fragment_atom = self.top_molecules[ligand][1]
                self.frag_core_atom = self.top_molecules[ligand][0]
                w.write(self.ligand_fragment)

    def get_top_molecule(self, heap, tmpdirname):
        """
        For each molecules within the top 100 most similar molecules from analyzed database:
            1. Assess if there is a substructure match between the initial seed compound and the external dataset molecule.
            2. Compute number of R-groups bound to the core.
            3. Extract PDB atom names from the linker atoms (the ones that form the bound between the core and the fragment).

        Parameters
        ----------
        heap: Heap with most similar molecules to input ligand. Maximum length is 100.
        tmpdirname: Temporal directory for temporal file storage.

        Returns
        -------
        output: Dictionary with RDkit molecule objects as keys, and PDB atom names of linker atoms as values.

        """
        print("Getting top molecules")
        no_substructure=0
        multiple_r_groups=0
        frag_core_atom_yaml = ""
        fragment_atom_yaml = ""
        output = {}
        for group in heap:
            molecule = Chem.MolFromPDBFile(os.path.join(tmpdirname, group[1] + ".pdb"), removeHs=True)
            if molecule:
                path = os.path.join(tmpdirname, group[1] + ".pdb")
                path_fragment = os.path.join(tmpdirname, group[1] + "_fragment.pdb")
                path_core = os.path.join(tmpdirname + "_core.pdb")

                mol = Chem.AddHs(molecule, addCoords=True)
                Chem.MolToPDBFile(mol, path)
                mol_with_hs = Chem.MolFromPDBFile(path, removeHs=False)

                mol_no_Hs = molecule
                mol_subs = molecule

                core_no_hs = self.init_mol

                self.substructure_results[group[1]] = mol_no_Hs.GetSubstructMatches(core_no_hs)
                # If there is a substructure match, convert to graph and count number of connected components
                if self.substructure_results[group[1]]:
                    bonds_no_Hs = [(m.GetBeginAtomIdx(), m.GetEndAtomIdx()) for m in mol_no_Hs.GetBonds()]
                    bonds = [(m.GetBeginAtomIdx(), m.GetEndAtomIdx()) for m in mol_with_hs.GetBonds()]

                    # Find number of connected components or each substructure match
                    h_len = []
                    for k in self.substructure_results[group[1]]:
                        H = self.generate_graph(bonds_no_Hs)
                        h_len.append(self.number_connected_components(H, k))

                    if 1 in h_len:
                        mw = Chem.RWMol(mol_subs)
                        mw2 = Chem.RWMol(mol_subs)

                        atoms_core =[]
                        atoms_fragment=[]
                        for atom in mw.GetAtoms():
                            if atom.GetIdx() not in self.substructure_results[group[1]][h_len.index(1)]:
                                atom.SetAtomicNum(0)
                                atoms_fragment.append(atom.GetIdx())

                        for atom in mw2.GetAtoms():
                            if atom.GetIdx() in self.substructure_results[group[1]][h_len.index(1)]:
                                atom.SetAtomicNum(0)
                                atoms_core.append(atom.GetIdx())

                        mw = Chem.DeleteSubstructs(mw, Chem.MolFromSmarts('[#0]'))
                        mw2 = Chem.DeleteSubstructs(mw2, Chem.MolFromSmarts('[#0]'))
                        fragment_no_hs = mw2
                        fragment_with_hs_nopdbinfo = Chem.AddHs(fragment_no_hs, addCoords=True)
                        Chem.MolToPDBFile(fragment_with_hs_nopdbinfo, path_fragment)
                        fragment_with_hs = Chem.MolFromPDBFile(path_fragment, removeHs=False)

                        core_with_hs_nopdbinfo = Chem.AddHs(mw, addCoords=True)
                        Chem.MolToPDBFile(core_with_hs_nopdbinfo, path_core)
                        core_with_hs = Chem.MolFromPDBFile(path_core, removeHs=False)

                        G = self.generate_graph(bonds)
                        counter = 0
                        for i, k in G.edges():
                            if (i in atoms_core and k in atoms_fragment) or (i in atoms_fragment and k in atoms_core):
                                counter += 1
                                if (i in atoms_core and k in atoms_fragment):
                                    core_atom_node = i
                                    fragment_atom_node = k
                                else:
                                    core_atom_node = k
                                    fragment_atom_node = i

                        _ = AllChem.GenerateDepictionMatching2DStructure(core_with_hs, self.init_mol_hs)
                        for atom in core_with_hs.GetAtoms():
                            if mol_with_hs.GetAtomWithIdx(core_atom_node).GetPDBResidueInfo().GetName().strip() == atom.GetPDBResidueInfo().GetName().strip():
                                coords_atom_mw = core_with_hs.GetConformers()[0].GetPositions()[atom.GetIdx()]

                        for atom_init in self.init_mol_hs.GetAtoms():
                            comparison = self.init_mol_hs.GetConformers()[0].GetPositions()[atom_init.GetIdx()][:2] == np.array(coords_atom_mw)[:2]
                            if comparison.all():
                                frag_core_atom_yaml = "%s" % (atom_init.GetPDBResidueInfo().GetName().strip())

                        for bond_frag in fragment_with_hs.GetBonds():
                            if bond_frag.GetBeginAtom().GetPDBResidueInfo().GetName().strip() == mol_with_hs.GetAtomWithIdx(fragment_atom_node).GetPDBResidueInfo().GetName().strip() and "H" in bond_frag.GetEndAtom().GetPDBResidueInfo().GetName().strip():
                                fragment_atom_yaml = "%s" % (bond_frag.GetEndAtom().GetPDBResidueInfo().GetName().strip())

                            elif bond_frag.GetEndAtom().GetPDBResidueInfo().GetName().strip() == mol_with_hs.GetAtomWithIdx(fragment_atom_node).GetPDBResidueInfo().GetName().strip() and "H" in bond_frag.GetBeginAtom().GetPDBResidueInfo().GetName().strip():
                                fragment_atom_yaml = "%s" % (bond_frag.GetBeginAtom().GetPDBResidueInfo().GetName().strip())

                        for bond in self.init_mol_hs.GetBonds():
                            if (bond.GetBeginAtom().GetPDBResidueInfo().GetName().strip() == frag_core_atom_yaml and "H" in bond.GetEndAtom().GetPDBResidueInfo().GetName().strip()):
                                frag_core_atom_yaml = "%s-%s" % (frag_core_atom_yaml, bond.GetEndAtom().GetPDBResidueInfo().GetName().strip())

                            if (bond.GetEndAtom().GetPDBResidueInfo().GetName().strip() == frag_core_atom_yaml and "H" in bond.GetBeginAtom().GetPDBResidueInfo().GetName().strip()):
                                frag_core_atom_yaml = "%s-%s" % (frag_core_atom_yaml, bond.GetBeginAtom().GetPDBResidueInfo().GetName().strip())

                        if counter <= 1:
                            output[group[1]] = [frag_core_atom_yaml, fragment_atom_yaml, fragment_with_hs]
                    else:
                        multiple_r_groups += 1
                else:
                    no_substructure+=1
        print("%s STRUCTURES DON'T HAVE THE INPUT LIGAND AS A SUBSTRUCTURE" % (no_substructure))
        print("%s STRUCTURES HAVE MORE THAN ONE R-GROUP" % (multiple_r_groups))

        return output

    def atoms_fragment(self, atoms_core, graph):
        H = graph.copy()
        H.remove_nodes_from(atoms_core)
        return list(H.nodes())

    def generate_graph(self, bonds):
        """
        For a given list of bonds, represented as relationships between atom indices, generate a graph of a molecule.

        Parameters
        ----------
        bonds: List of tuples. Each tuple consists of two atom indices.

        Returns
        -------
        G: Graph of the molecule.

        """
        G = nx.Graph()
        for i in bonds:
            G.add_edge(i[0], i[1])
        return G

    def number_connected_components(self, graph, substructure):
        for j in substructure:
            graph.remove_node(j)
        return nx.number_connected_components(graph)

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
        """
        If a yaml file with specified filters is provided by the user,
        assess if a molecule passes those filters.

        Returns
        -------
        bool: True if molecule fits filters, False otherwise.
        """
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
                    return True
        return False


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
    lib = Library(path, ligand, ligand_path, filters)
    return lib.top_molecules