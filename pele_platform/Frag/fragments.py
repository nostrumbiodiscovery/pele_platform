class Fragment:

    def __init__(self, file, atom_fragment, hydrogen_fragment, atom_core, hydrogen_core, no_hydrogen):
        from rdkit import Chem
        self.file = file
        self.atom_fragment = atom_fragment
        self.hydrogen_fragment = hydrogen_fragment
        self.atom_core = atom_core
        self.hydrogen_core = hydrogen_core
        self.mol = Chem.MolFromPDBFile(self.file, removeHs=False)
        self.no_hydrogen = no_hydrogen

    def get_inputfile_line(self):
        if self.no_hydrogen:
            self.line = "{} {} {}".format(self.file, 
                self.atom_core.get_name(),
                self.atom_fragment.get_name())
        else:
            self.line = "{} {}-{} {}-{}".format(self.file, 
                self.atom_core.get_name(), self.hydrogen_core.get_name(),
                self.atom_fragment.get_name(), self.hydrogen_fragment.get_name())
        return self.line

    def sanitize_file(self):
        with open(self.file, "r") as fin:
            lines = fin.readlines()
        with open(self.file, "w") as fout:
            new_lines = [ line[:17] + "FRG L" + line[22:] if line.startswith("HETATM") else line for line in lines]
            fout.write("".join(new_lines))
