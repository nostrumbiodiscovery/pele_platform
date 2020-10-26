from dataclasses import dataclass
from rdkit import Chem
import pele_platform.Errors.custom_errors as ce

@dataclass
class SmilesConstraints:

    input_pdb: str
    constrain_smiles: str
    resname: str
    chain: str
    spring_constant: float = 50.0

    def run(self):
        self._convert_to_smarts()
        self._extract_ligand()
        self._substructure_search()
        if self.substructures:
            self._build_constraints()
            return self.constraints
        else:
            raise ce.SubstructureError("Could not find {} substructure in the ligand {}.".format(self.constrain_smiles, Chem.MolToSmiles(self.ligand)))
    
    def _convert_to_smarts(self):
        self.pattern = Chem.MolFromSmiles(self.constrain_smiles)
        Chem.SanitizeMol(self.pattern)
        Chem.MolToSmarts(self.pattern)

    def _extract_ligand(self):
        complex = Chem.MolFromPDBFile(self.input_pdb)
        self.ligand = Chem.rdmolops.SplitMolByPDBResidues(complex)[self.resname]
    
    def _substructure_search(self):
        self.substructures = self.ligand.GetSubstructMatches(self.pattern)
        if not self.substructures:
            # try to convert molecule to SMARTS, manually remove aromatic bonds and load again as substructure pattern
            self.pattern = Chem.MolToSmarts(self.pattern).replace(":", "-")
            self.pattern = Chem.MolFromSmarts(self.pattern)
            self.substructures = self.ligand.GetSubstructMatches(self.pattern)

    def _build_constraints(self):
        self.constraints = []
        if len(self.substructures) > 1:
            raise ce.SubstructureError("More than one substructure found in your ligand. Make sure SMILES constrain pattern is not ambiguous!")
        else:
            for m in self.substructures[0]:
                atom = self.ligand.GetAtomWithIdx(m).GetMonomerInfo()
                template = '{{ "type": "constrainAtomToPosition", "springConstant": {}, "equilibriumDistance": 0.0, "constrainThisAtom": "{}:{}:{}" }},'
                self.constraints.append(template.format(self.spring_constant, self.chain, atom.GetResidueNumber(), atom.GetName().strip()))
