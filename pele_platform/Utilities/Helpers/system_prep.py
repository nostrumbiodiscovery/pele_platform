import os
import sys
import subprocess
import argparse

sys.path.append(
    os.path.dirname(
        os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    )
)
import pele_platform.constants.constants as cs


class SystemBuilder(object):
    def __init__(self, receptor, ligand, residue, pele_dir, inputs_dir):
        self.receptor = receptor
        self.ligand = ligand
        self.residue = residue
        self.pele_dir = pele_dir
        self.system = self.receptor
        self.inputs_dir = inputs_dir

    @classmethod
    def build_system(cls, receptor, ligand, residue, pele_dir, output=False, inputs_dir=None):
        SPYTHON = os.path.join(cs.SCHRODINGER, "utilities/python")
        if not os.path.exists(SPYTHON):
            SPYTHON = os.path.join(cs.SCHRODINGER, "run")
        system = cls(receptor, ligand, residue, pele_dir, inputs_dir)
        system.receptor, system.lig_ref = system.retrieve_receptor(output=output)
        system.lig = ligand if ligand else "{}.mae".format(residue)
        my_env = os.environ.copy()
        my_env["PYTHONPATH"] = ""
        my_env["SCHRODINGER_PYTHONPATH"] = os.path.join(
            cs.SCHRODINGER, "internal/lib/python2.7/site-packages/"
        )
        my_env["SCHRODINGER"] = cs.SCHRODINGER
        subprocess.call(
            "{} {} {} {}".format(SPYTHON, __file__, system.lig_ref, pele_dir).split(),
            env=my_env,
        )
        system.residue = residue
        return system

    def build_complex(self):
        """
        From the receptor and ligand in pdb build
        another pdb with the whole complex
        """
        complex_content = []

        name = os.path.basename(os.path.splitext(self.receptor)[0])
        self.complex = os.path.join(self.inputs_dir, "{}_complex.pdb".format(name))

        try:
            with open(self.receptor, "r") as pdb_file:
                receptor_text = [
                    line
                    for line in pdb_file
                    if line.startswith("ATOM")
                    or line.startswith("HETATM")
                    or line.startswith("TER")
                ]
            with open(self.lig_ref, "r") as pdb_file:
                ligand_text = [line for line in pdb_file if line.startswith("HETATM")]
        except IOError:
            raise IOError(
                "Receptor or ligand not found. Check your ligand residue/chain name or input files"
            )

        if not receptor_text or not ligand_text:
            raise ValueError(
                "The ligand_pdb was not properly created. Check your MAE file."
            )

        complex_content.extend(receptor_text + ["TER\n"] + ligand_text + ["END"])

        with open(self.complex, "w") as fout:
            fout.write("".join(complex_content))

        return self.complex

    def convert_mae(self):
        """
        Desciption: From each structure retrieve
        a .mae file of the ligand in the receptor.
        Output:
             structure_mae: ligand
             res = residue
        """

        for structure in st.StructureReader(self.lig_ref):
            for residue in structure.residue:
                res = residue.pdbres.strip()
            str_name = "{}".format(res)
            try:
                structure.write(str_name + ".mae")
            except ValueError:
                str_name = "{}".format(res)
            finally:
                structure.write(str_name + ".mae")
                structure_mae = "{}.mae".format(str_name)
        return structure_mae, res

    def retrieve_receptor(self, output=False):
        """
        This function returns receptor of the complex of interest.
        :param complex: system format pdb
        :output: receptor text
        """
        ligand = output if output else os.path.join(self.pele_dir, "input", "ligand.pdb")
        receptor = os.path.join(self.inputs_dir, "receptor.pdb")

        with open(self.receptor, "r") as pdb_file:
            receptor_text = [
                line
                for line in pdb_file
                if line.startswith("ATOM")
                or "TER" in line
                or (line.startswith("HETATM") and line[17:20].strip() != self.residue)
            ]
        with open(self.receptor, "r") as pdb_file:
            residue = None
            ligand_text = []

            # Filter out all description, CONECT and ANISOU lines
            lines = [line for line in pdb_file.readlines() if line.startswith("ATOM") or line.startswith("HETATM")]

            for line in lines:
                if line[17:20].strip() == self.residue:
                    if residue is None:
                        residue = line[22:26]
                        chain = line[20:22]
                        ligand_text.append(line)
                    else:
                        if residue == line[22:26] and chain == line[20:22]:
                            ligand_text.append(line)

        if not receptor_text or not ligand_text:
            raise ValueError(
                "Something went wrong when extracting the ligand. Check residue&Chain on input"
            )
        with open(ligand, "w") as fout:
            fout.write("".join(ligand_text + ["END"]))
        with open(receptor, "w") as fout:
            fout.write("".join(receptor_text + ["END"]))

        return "".join(receptor_text), ligand


def convert_pdb(mae_file, output_dir):
    from schrodinger import structure as st

    for structure in st.StructureReader(mae_file):
        structure.write(output_dir)


def convert_mae(pdb, output_dir="."):
    """
    Desciption: From each structure retrieve
    a .mae file of the ligand in the receptor.
    Output:
         structure_mae: ligand
         res = residue
    """
    from schrodinger import structure as st

    for structure in st.StructureReader(pdb):
        for residue in structure.residue:
            res = residue.pdbres.strip()
            str_name = os.path.abspath(os.path.join(output_dir, "{}".format(res)))
            try:
                structure.write(str_name + ".mae")
            except ValueError:
                str_name = "{}".format(res)
            finally:
                structure.write(str_name + ".mae")
                structure_mae = "{}.mae".format(str_name)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("ligand", help="ligand input file to convert")
    parser.add_argument(
        "output_dir", help="output directory to dump the converted file"
    )
    parser.add_argument(
        "--mae",
        action="store_true",
        help="Whether to convert to mae (--mae) or pdb (not --mae)",
    )
    args = parser.parse_args()
    return args.ligand, args.output_dir, args.mae


if __name__ == "__main__":
    input_file, output_dir, ligand_mae = parse_args()
    if ligand_mae:
        convert_pdb(input_file, output_dir)
    else:
        convert_mae(input_file, output_dir)
