from dataclasses import dataclass
import re
import pele_platform.Errors.custom_errors as ce
import pele_platform.Utilities.Helpers.helpers as hp


@dataclass
class SmilesConstraints:
    input_pdb: str
    constrain_core: str
    resname: str
    chain: str
    spring_constant: float = 50.0

    def run(self):
        from rdkit import Chem
        self.smarts = self.convert_to_smarts(self.constrain_core)
        self.ligand = self.extract_ligand(self.input_pdb, self.resname)
        self.substructures = self.get_matches(self.smarts, self.ligand)
        if self.substructures:
            self.constraints = self.build_constraints(self.substructures, self.ligand, self.spring_constant, self.chain)
            return self.constraints
        else:
            raise ce.SubstructureError(
                "Could not recognise the substructure specified in 'constrain_core'. This might be due to differences "
                "in ionisation states. Check, if the charges are identical and adjust your pattern, if necessary.\n"
                "Core SMARTS: {}\nLigand SMARTS: {}".format(
                    self.smarts, Chem.MolToSmarts(self.ligand)))

    def convert_to_smarts(self, core):
        from rdkit import Chem
        if "C" in core or "c" in core:  # if the pattern is SMILES
            pattern = Chem.MolFromSmiles(core)
            smarts = Chem.MolToSmarts(pattern)
        else:  # if the pattern is SMARTS
            smarts = core
        return smarts

    def extract_ligand(self, pdb, resname):
        from rdkit import Chem
        complex = Chem.MolFromPDBFile(pdb)
        if complex is not None:
            ligand = Chem.rdmolops.SplitMolByPDBResidues(complex)[resname]
        else:
            ligand = self._backup_ligand_extraction(pdb, resname)
        return ligand

    def get_matches(self, smarts, ligand):
        substructures = self._substructure_search(smarts, ligand)
        if not substructures:
            substructures = self._substructure_search(smarts, ligand, ":", "-")  # removing aromatic bonds
            if not substructures:
                substructures = self._substructure_search(smarts, ligand, "=", "-")  # removing double bonds
                if not substructures:
                    substructures = self._substructure_search(smarts, ligand, regex_remove=True)  # remove hydrogens with regex
                    if not substructures:
                        substructures = self._substructure_search(smarts, ligand,
                            rdkit_remove=True)  # remove hydrogens by iterating through all atoms
        return substructures

    def _substructure_search(self, smarts, ligand, old="", new="", regex_remove=False, rdkit_remove=False):
        from rdkit import Chem
        self.smarts = smarts
        if regex_remove:
            self.smarts = re.sub(r"H\d?", "", self.smarts)
        else:
            self.smarts = self.smarts.replace(old, new)
        self.pattern = Chem.MolFromSmarts(self.smarts)

        if rdkit_remove:
            idx_to_remove = []
            for atom in self.pattern.GetAtoms():
                if atom.GetSymbol() == "H":
                    idx_to_remove.append(atom.GetIdx())
            for i in idx_to_remove:
                self.pattern.RemoveAtom(i)
        return ligand.GetSubstructMatches(self.pattern)

    def build_constraints(self, substructures, ligand, spring_constant, chain):
        from rdkit import Chem
        constraints = []
        if len(substructures) > 1:
            raise ce.SubstructureError(
                "More than one substructure found in your ligand. Make sure SMILES constrain pattern is not ambiguous!")
        else:
            for m in substructures[0]:
                atom = ligand.GetAtomWithIdx(m).GetMonomerInfo()
                template = '{{ "type": "constrainAtomToPosition", "springConstant": {}, "equilibriumDistance": 0.0, "constrainThisAtom": "{}:{}:{}" }},'
                constraints.append(template.format(spring_constant, chain, atom.GetResidueNumber(),
                                                        atom.GetName().replace(" ", "_")))
        return constraints

    def _backup_ligand_extraction(self, pdb, resname):
        """
        Extracts ligand lines from PDB file and loads them directly. Sometimes rdkit refuses to read in PDB files
        because of random valence errors, hence the need for a backup if extract_ligand fails.
        """
        from rdkit import Chem
        ligand_lines = []

        with open(pdb, "r") as f:
            lines = f.readlines()
            for line in lines:
                if (line.startswith("ATOM") or line.startswith("HETATM")) and line[17:20].strip() == resname:
                    ligand_lines.append(line)
        ligand_block = "".join(ligand_lines)
        ligand = Chem.MolFromPDBBlock(ligand_block)
        return ligand
