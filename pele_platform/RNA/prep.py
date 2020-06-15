import os
import argparse



class PDBFix():

    def __init__(self, pdb: str, output_folder=":"):
        self.pdb = pdb
        self.lines = self.get_lines()
        self.pdb_out = self.pdb.rsplit(".", 1)[0] + "_proc.pdb"
        self.output_folder = output_folder

    def get_lines(self):
        with open(self.pdb, "r") as f:
            return f.readlines()

    def fix_nucleotides(self):
        new_lines = []
        for line in self.lines:
            try:
                atom_name = line[12:16].strip()
            except IndexError:
                new_lines.append(line)
            if atom_name == "H5'1":
                atom_name = "H5''"
            elif atom_name == "H5'2":
                atom_name = " H5'"
            elif atom_name == "H6 1":
                atom_name = " H61"
            elif atom_name == "H6 2":
                atom_name = " H62"
            elif atom_name == "HN1":
                atom_name = " H1 "
            elif atom_name == "HN3":
                atom_name = " H3 "
            elif atom_name == "H4 1":
                atom_name = " H41"
            elif atom_name == "H4 2":
                atom_name = " H42"
            elif atom_name == "H2 1":
                atom_name = " H21"
            elif atom_name == "H2 2":
                atom_name = " H22"
            else:
                atom_name = line[12:16]
            new_line = line[:12] + atom_name + line[16:]
            new_lines.append(new_line)
        self.lines = new_lines

    def change_K_to_na(self):
        self.lines = [line.replace(" K ", "NA ") if line.startswith("HETATM") else line for line in self.lines]


    def remove_two_first(self):
        residues_to_remove = [line[21:26] for i, line in enumerate(self.lines) if line.startswith("ATOM") and not (self.lines[i+1].startswith("ATOM") and self.lines[i-1].startswith("ATOM"))]
        new_lines = []
        for line in self.lines:
            if line.startswith("ATOM"):
                residue = line[21:26]
            else:
                new_lines.append(line)
                continue
            if residue not in residues_to_remove:
                new_lines.append(line)
        self.lines = new_lines

    def write(self):
        path = os.path.join(self.output_folder, self.pdb_out)
        with open(path, "w") as f:
            f.write("".join(self.lines))
        return path

def fix_rna_pdb(pdb):
    env.logger("Fixing pdb for RNA")
    pdb = PDBFix(pdb)
    pdb.fix_nucleotides()
    pdb.remove_two_first()
    pdb.change_K_to_na()
    output = pdb.write()
    return output

def parse_args(parser):
    parser.add_argument('pdb', type=str, help="Complex pdb with ligand-RNA")
    return parser


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='launch DNA simulation')
    parser = parse_args(parser)
    args = parser.parse_args()
    pdb_fixed = main(args.pdb)
            
               
