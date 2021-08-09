import os

from pele_platform.Adaptive.parametrizer import Parametrizer
from peleffy import topology


class LigandConformations:
    def __init__(self, path, system, resname, forcefield, pele_dir=None):
        """
        Initializes the LigandConformations class.

        Parameters
        ----------
        path : str
            Path to the folder with ligand conformations in PDB format
        system : str
            Path to system PDB file.
        forcefield : str
            Forcefield to parametrize the ligand, as set by the user in YAML.
        pele_dir : str
            Path to save conformations library.
        """
        self.path = path
        self.resname = resname
        self.ligand = self._extract_ligand(system)
        self.forcefield = Parametrizer._retrieve_forcefield(forcefield)
        self.dir = pele_dir if pele_dir is not None else os.getcwd()
        self.output_file = os.path.join(
            self.dir, "DataLocal", "Conformations", f"{resname}.conformation"
        )

    def generate(self):
        """
        Generates conformations library.

        Returns
        -------
            Path to LIG.conformation file.
        """
        self._parametrize_ligand()
        output = self._get_conformations()
        return output

    def _parametrize_ligand(self):
        """
        Parametrizes the ligand based on user-defined forcefield and generates a peleffy Topology object.
        """
        self.ligand_parameters = self.forcefield.parameterize(self.ligand)
        self.topology = topology.Topology(self.ligand, self.ligand_parameters)

    def _get_conformations(self):
        """
        Calculates conformations of the ligand and saves them to DataLocal/conformations folder.
        """
        bce = topology.BCEConformations(self.topology, self.path, from_bce=False)
        bce.calculate()

        dir_tree = os.path.dirname(self.output_file)
        if not os.path.exists(dir_tree):
            os.makedirs(dir_tree)

        bce.save(self.output_file)
        return self.output_file

    def _extract_ligand(self, system):
        """
        Extracts ligand lines and its CONECT lines from system PDB.

        Parameters
        ----------
        system : str
            Path to system PDB file.

        Returns
        -------
            Peleffy Molecule object with the ligand.
        """
        temp_file = "temp_ligand.pdb"

        with open(system, "r") as file:
            lines = file.readlines()

        pdb_lines = [
            line
            for line in lines
            if line.startswith("ATOM") or line.startswith("HETATM")
        ]
        conect_lines = [line for line in lines if line.startswith("CONECT")]

        ligand_lines = list()
        for line in pdb_lines:
            if line[17:20].strip() == self.resname:
                ligand_lines.append(line)

        with open(temp_file, "w+") as temp:
            for line in ligand_lines:
                temp.write(line)
            for line in conect_lines:
                temp.write(line)

        mol = topology.Molecule(temp_file, allow_undefined_stereo=True)
        os.remove(temp_file)
        return mol
