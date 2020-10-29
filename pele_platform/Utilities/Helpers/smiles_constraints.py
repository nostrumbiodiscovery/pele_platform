from dataclasses import dataclass
import re
import pele_platform.Errors.custom_errors as ce

@dataclass
class SmilesConstraints:

    input_pdb: str
    constrain_core: str
    resname: str
    chain: str
    spring_constant: float = 50.0

    def run(self):
        from rdkit import Chem
        self._convert_to_smarts()
        self._extract_ligand()
        self._get_matches()
        if self.substructures:
            self._build_constraints()
            return self.constraints
        else:
            raise ce.SubstructureError("Could not recognise the substructure specified in 'constrain_core'. This might be due to differences in ionisation states. Check, if the charges are identical and adjust your pattern, if necessary.\nCore SMARTS: {}\nLigand SMARTS: {}".format(self.smarts, Chem.MolToSmarts(self.ligand)))
    
    def _convert_to_smarts(self):
        from rdkit import Chem
        if ("C" or "c") in self.constrain_core:  # if the pattern is SMILES
            self.pattern = Chem.MolFromSmiles(self.constrain_core)
            self.smarts = Chem.MolToSmarts(self.pattern)
        else:  # if the pattern is SMARTS
            self.pattern = Chem.MolFromSmarts(self.constrain_core)
            self.smarts = self.constrain_core

    def _extract_ligand(self):
        from rdkit import Chem
        complex = Chem.MolFromPDBFile(self.input_pdb)
        if complex is not None:
            self.ligand = Chem.rdmolops.SplitMolByPDBResidues(complex)[self.resname]
        else:
            self._backup_ligand_extraction()

    def _get_matches(self):
        self.substructures = self._substructure_search()
        if not self.substructures:
            self.substructures = self._substructure_search(":", "-")  # removing aromatic bonds
            if not self.substructures:
                self.substructures = self._substructure_search("=", "-")  # removing double bonds
                if not self.substructures:
                    self.substructures = self._substructure_search(regex_remove=True)  # remove hydrogens with regex
                    if not self.substructures:
                        self.substructures = self._substructure_search(rdkit_remove=True)  # remove hydrogens by iterating through all atoms

    def _substructure_search(self, old="", new="", regex_remove=False, rdkit_remove=False):
        from rdkit import Chem
        if regex_remove:
            self.smarts = re.sub(r"H\d?","",self.smarts)
        else:
            self.smarts = self.smarts.replace(old, new)
        print("Trying SMARTS", self.smarts)        
        self.pattern = Chem.MolFromSmarts(self.smarts) 
        
        if rdkit_remove:
            print("RDkit")
            idx_to_remove = []
            for atom in self.pattern.GetAtoms():
                if atom.GetSymbol() == "H":
                    idx_to_remove.append(atom.GetIdx())
            for i in idx_to_remove:
                self.pattern.RemoveAtom(i)
            print(len(self.pattern.GetAtoms()))

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
        from rdkit import Chem
        ligand_lines = []

        with open(self.input_pdb, "r") as f:
            lines = f.readlines()
            for line in lines:
                if (line.startswith("ATOM") or line.startswith("HETATM")) and line[17:20].strip() == self.resname:
                    ligand_lines.append(line)

        ligand_block = "".join(ligand_lines)
        self.ligand = Chem.MolFromPDBBlock(ligand_block)
