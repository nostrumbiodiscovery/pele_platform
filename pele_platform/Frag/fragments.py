




class Fragment():


    def __init__(self, file, atom_fragment_idx, hydrogen_added_idx, atom_core_name, hydrogen_core_name, no_hydrogen):
        from rdkit import Chem
        import rdkit.Chem.rdchem as rc
        self.file = file
        self.atom_fragment_idx = atom_fragment_idx
        self.hydrogen_added_idx = hydrogen_added_idx
        self.atom_core_name = atom_core_name.strip()
        self.hydrogen_core_name = hydrogen_core_name.strip()
        self.mol = Chem.MolFromPDBFile(self.file, removeHs=False)
        self.atom_fragment_name = self.mol.GetAtomWithIdx(self.atom_fragment_idx).GetMonomerInfo().GetName().strip()
        self.atom_hydrogen_name = self.mol.GetAtomWithIdx(self.hydrogen_added_idx).GetMonomerInfo().GetName().strip()
        self.no_hydrogen = no_hydrogen


    def get_inputfile_line(self):
        if self.no_hydrogen:
            self.line = "{} {} {}".format(self.file, 
                self.atom_core_name,
                self.atom_fragment_name)
        else:
            self.line = "{} {}-{} {}-{}".format(self.file, 
                self.atom_core_name, self.hydrogen_core_name,
                self.atom_fragment_name, self.atom_hydrogen_name)
        return self.line

    def santize_file(self):
        with open(self.file, "r") as fin:
            lines = fin.readlines()
        with open(self.file, "w") as fout:
            new_lines = [ line[:17] + "FRG L" + line[22:] if line.startswith("HETATM") else line for line in lines]
            fout.write("".join(new_lines))
