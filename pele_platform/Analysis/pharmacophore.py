import pele_platform.Utilities.Helpers.helpers as hp


class PharmacophoreAnalysis():

    def __init__(self):
        pass

    def run(self):
        pass
        



class RdKitPharmacoPhore():
    """
    Base class with utilities to retrieve 
    pharmacophore-atom relationship
    """

    def __init__(self, pdb):
        if hp.is_rdkit():
            from rdkit import RDConfig
            import os
            from rdkit.Chem import ChemicalFeatures
        self.pdb = pdb
        fdefName = os.path.join(RDConfig.RDDataDir,'BaseFeatures.fdef')
        self.factory = ChemicalFeatures.BuildFeatureFactory(fdefName)
        self.structure = None
        self.feats = None

    def get_pharmacophore_atom_realtionship(self):
        structure = self.load_mol()
        self.hds = self.get_hd()
        self.has = self.get_ha()
        self.aromatic = self.get_aromatic()
        self.hydro = self.get_hydro()
        print(len(self.hds), len(self.has), len(self.aromatic), len(self.hydro))
        self.draw_pharmacophore()

    def load_mol(self):
        from rdkit import Chem
        from rdkit.Chem import rdmolfiles
        if self.pdb[-4:] == ".pdb":
            self.structure =  rdmolfiles.MolFromPDBFile(self.pdb, removeHs=False)
        elif self.pdb[-4:] == ".sdf":
            self.structure = next(Chem.SDMolSupplier(self.pdb, removeHs=False))
        return self.structure

    def get_hd(self):
        from rdkit import Chem
        hbonds = []
        for atom in self.structure.GetAtoms():
            if atom.GetSymbol() in ["N", "O"] and [a for a in atom.GetNeighbors() if a.GetSymbol() == "H"]:
                hbonds.append(atom)
        return hbonds

    def get_ha(self):
        from rdkit import Chem
        hbonds = []
        for atom in self.structure.GetAtoms():
            if atom.GetSymbol() in ["N", "O"] and not [a for a in atom.GetNeighbors() if a.GetSymbol() == "H"] and atom.GetImplicitValence() != 0:
                hbonds.append(atom)
        return hbonds

    def get_aromatic(self):
        from rdkit import Chem
        aroms = []
        for atom in self.structure.GetAtoms():
            if atom.GetSymbol() == "C" and atom.IsInRing() and atom.GetTotalNumHs() ==1:
                aroms.append(atom)
        return aroms

    def get_hydro(self):
        hydr = []
        for atom in self.structure.GetAtoms():
            if atom.GetSymbol() in ["F", "Cl","Br"]:
                hydr.append(atom)
            elif atom.GetSymbol() == "C" and len([a for a in atom.GetNeighbors() if a.GetSymbol() == "H"]) > 1:
                hydr.append(atom)
        return hydr

    def draw_pharmacophore(self):
        lines = []
        pharm = {"O": self.hds, "H": self.has, "C": self.aromatic, "F": self.hydro}
        coords = self.structure.GetConformers()[0].GetPositions()
        for element, atoms in pharm.items():
            for atom in atoms:
                idx = atom.GetIdx()
                lines.append(format_line_pdb(coords[idx], element, bfact=1))
        with open(f"{self.pdb[:-4]}_pharmacophore.pdb", "w") as f:
            f.write("".join(lines))
                


def format_line_pdb(coords, atomname, bfact, atomnum = "1", resname="UNK", chain="A", resnum="1", occ=1.00):
    x, y, z = coords
    atomstr = "ATOM"
    element = atomname[0]
    line = f"{atomstr:5}{atomnum:>5} {atomname:4} {resname:3} {chain}{resnum:>4}    {x:>8.3f}{y:>8.3f}{z:>8.3f}{occ:6.2f}{bfact:6.2f}{element:>12}\n"
    return line


if __name__ == "__main__":
    pharm = RdKitPharmacoPhore("ligand.pdb")
    pharm.get_pharmacophore_atom_realtionship()

        
        
        


     
        
