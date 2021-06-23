class Atom:

    def __init__(self, molecule, idx):
        self.molecule = molecule
        self.idx = idx

    def get_name(self):
        self.name = self.molecule.GetAtomWithIdx(self.idx).GetMonomerInfo().GetName().strip()
        return self.name
