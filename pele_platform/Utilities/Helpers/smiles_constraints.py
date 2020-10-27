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
        if complex is not None:
            self.ligand = Chem.rdmolops.SplitMolByPDBResidues(complex)[self.resname]
        else:
            self._backup_ligand_extraction()

    def _substructure_search(self):
        print("Ligand SMARTS", Chem.MolToSmarts(self.ligand))
        self.substructures = self.ligand.GetSubstructMatches(self.pattern)
        print("Original pattern", Chem.MolToSmarts(self.pattern))
        if not self.substructures:
            self.substructures = self._change_smarts(":", "-")  # removing aromatic bonds
            if not self.substructures:
                self.substructures = self._change_smarts("=", "-")  # removing double bonds
                if not self.substructures:
                    self.substructures = self._change_smarts("@", "")  # removing stereochemistry

    def _change_smarts(self, old, new):
        self.pattern = Chem.MolToSmarts(self.pattern).replace(old, new)
        print("Trying pattern", self.pattern)
        self.pattern = Chem.MolFromSmarts(self.pattern)
        return self.ligand.GetSubstructMatches(self.pattern)

    def _build_constraints(self):
        self.constraints = []
        if len(self.substructures) > 1:
            raise ce.SubstructureError("More than one substructure found in your ligand. Make sure SMILES constrain pattern is not ambiguous!")
        else:
            for m in self.substructures[0]:
                atom = self.ligand.GetAtomWithIdx(m).GetMonomerInfo()
                template = '{{ "type": "constrainAtomToPosition", "springConstant": {}, "equilibriumDistance": 0.0, "constrainThisAtom": "{}:{}:{}" }},'
                self.constraints.append(template.format(self.spring_constant, self.chain, atom.GetResidueNumber(), atom.GetName().replace(" ", "_")))

    def _backup_ligand_extraction(self):
        ligand_lines = []

        with open(self.input_pdb, "r") as f:
            lines = f.readlines()
            for line in lines:
                if (line.startswith("ATOM") or line.startswith("HETATM")) and line[17:20].strip() == self.resname:
                    ligand_lines.append(line)

        ligand_block = "".join(ligand_lines)
        self.ligand = Chem.MolFromPDBBlock(ligand_block)
