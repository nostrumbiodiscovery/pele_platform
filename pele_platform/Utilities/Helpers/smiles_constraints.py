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
    hp.is_rdkit()

    def run(self):
        self.smarts = self._convert_to_smarts()
        self.ligand = self._extract_ligand()
        self.substructures = self._get_matches()
        if self.substructures:
            self.constraints = self._build_constraints()
            return self.constraints
        else:
            raise ce.SubstructureError(
                "Could not recognise the substructure specified in 'constrain_core'. This might be due to differences "
                "in ionisation states. Check, if the charges are identical and adjust your pattern, if necessary.\n"
                "Core SMARTS: {}\nLigand SMARTS: {}".format(
                    self.smarts, Chem.MolToSmarts(self.ligand)))

    def _convert_to_smarts(self):
        if "C" in self.constrain_core or "c" in self.constrain_core:  # if the pattern is SMILES
            pattern = Chem.MolFromSmiles(self.constrain_core)
            smarts = Chem.MolToSmarts(pattern)
        else:  # if the pattern is SMARTS
            smarts = self.constrain_core
        return smarts

    def _extract_ligand(self):
        complex = Chem.MolFromPDBFile(self.input_pdb)
        if complex is not None:
            ligand = Chem.rdmolops.SplitMolByPDBResidues(complex)[self.resname]
        else:
            ligand = self._backup_ligand_extraction()
        return ligand

    def _get_matches(self):
        substructures = self._substructure_search()
        if not substructures:
            substructures = self._substructure_search(":", "-")  # removing aromatic bonds
            if not substructures:
                substructures = self._substructure_search("=", "-")  # removing double bonds
                if not substructures:
                    substructures = self._substructure_search(regex_remove=True)  # remove hydrogens with regex
                    if not substructures:
                        substructures = self._substructure_search(
                            rdkit_remove=True)  # remove hydrogens by iterating through all atoms
        return substructures

    def _substructure_search(self, old="", new="", regex_remove=False, rdkit_remove=False):
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

        return self.ligand.GetSubstructMatches(self.pattern)

    def _build_constraints(self):
        constraints = []
        if len(self.substructures) > 1:
            raise ce.SubstructureError(
                "More than one substructure found in your ligand. Make sure SMILES constrain pattern is not ambiguous!")
        else:
            for m in self.substructures[0]:
                atom = self.ligand.GetAtomWithIdx(m).GetMonomerInfo()
                template = '{{ "type": "constrainAtomToPosition", "springConstant": {}, "equilibriumDistance": 0.0, "constrainThisAtom": "{}:{}:{}" }},'
                constraints.append(template.format(self.spring_constant, self.chain, atom.GetResidueNumber(),
                                                        atom.GetName().replace(" ", "_")))

        return constraints

    def _backup_ligand_extraction(self):
        """
        Extracts ligand lines from PDB file and loads them directly. Sometimes rdkit refuses to read in PDB files
        because of random valence errors, hence the need for a backup if _extract_ligand fails.
        """
        ligand_lines = []

        with open(self.input_pdb, "r") as f:
            lines = f.readlines()
            for line in lines:
                if (line.startswith("ATOM") or line.startswith("HETATM")) and line[17:20].strip() == self.resname:
                    ligand_lines.append(line)

        ligand_block = "".join(ligand_lines)
        ligand = Chem.MolFromPDBBlock(ligand_block)
        return ligand
