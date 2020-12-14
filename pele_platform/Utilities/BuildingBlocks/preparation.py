import os
from pele_platform.Utilities.Parameters import pele_env as pv


class Preparation:
    simulation_params: pv.EnviroBuilder

    def run(self, remove_water=False):

        to_remove = []

        # remove additional protein chains and all water molecules
        with open(self.simulation_params.system, "r") as file:
            lines = file.readlines()

            for line in lines:
                if ((line.startswith("ATOM") or line.startswith("HETATM")) and line[
                        21:22].strip() not in self.simulation_params.protein) or line.startswith(
                        "END") or line.startswith("TER") or line.startswith("CONECT"):
                    to_remove.append(line)
                if remove_water:
                    if (line.startswith("ATOM") or line.startswith("HETATM")) and line[17:20].strip() == "HOH":
                        to_remove.append(line)

            protein = [line for line in lines if line not in to_remove]

        # read in ligand file
        ligand = []
        with open(self.simulation_params.ligand_pdb, "r") as ligand_file:
            lines = ligand_file.readlines()
            for line in lines:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    ligand.append(line)

        self.simulation_params.system = os.path.abspath(
            os.path.basename(self.simulation_params.system).replace(".pdb", "_prep.pdb"))

        # join protein and ligand PDBs into new file
        with open(self.simulation_params.system, "w+") as file:
            for line in protein:
                file.write(line)
            file.write("\n")
            for line in ligand:
                file.write(line)

        return self.simulation_params
