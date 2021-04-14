from peleffy.utils.input import PDB
from peleffy.utils import OutputPathHandler
from peleffy.topology import RotamerLibrary, Topology
from peleffy import forcefield as ff
from peleffy.template import Impact

from pele_platform.constants import constants
from pele_platform.Errors import custom_errors


class LigandParametrization:
    """
    Base class to generate ligand parameters:
    1) Generates forcefield and rotamers from ligand.
    2) Copies user's external files.
    """

    TEMPLATE_FILE = "{}z"
    ROTAMER_FILE = "{}.rot.assign"

    # Mapping between args.forcefield and peleffy classes
    forcefields = {
        "opls2005": ff.OPLS2005ForceField(),
        "openff-1.3.0": ff.OpenForceField("openff_unconstrained-1.3.0.offxml"),
        "openff-1.2.1": ff.OpenForceField("openff_unconstrained-1.2.1.offxml"),
        "openff-1.2.0": ff.OpenForceField("openff_unconstrained-1.2.0.offxml"),
        "openff-1.1.1": ff.OpenForceField("openff_unconstrained-1.1.1.offxml"),
        "openff-1.1.0": ff.OpenForceField("openff_unconstrained-1.1.0.offxml"),
        "openff-1.0.1": ff.OpenForceField("openff_unconstrained-1.0.1.offxml"),
        "openff-1.0.0": ff.OpenForceField("openff_unconstrained-1.0.0.offxml"),
    }

    charge_parametrization_methods = ["am1bcc", "gasteiger", "OPLS"]

    def __init__(
        self,
        pdb_file,
        forcefield="OPLS2005",
        charge_parametrization_method="am1bcc",
        gridres=10,
        solvent="OBC",
        external_templates=None,
        external_rotamers=None,
        as_datalocal=True,
        pele_dir=None,
    ):
        """
        Initializes LigandParametrization to generate template and rotamer files.

        Parameters
        ------------
        pdb_file : str
            Path to the PDB file from which all HET groups will be extracted and parametrized.
        forcefield : str
            User-defined string to select forcefield for parametrization. Default = "OPLS2005".
        charge_parametrization_method : str
            User-defined string to select charge parametrization method. Default = "am1bcc".
        gridres : int
            Resolution of the rotamers when sampling. Default = 10 degrees.
        solvent : str
            Simulation solvent. Default = "OBC" if using OpenFF, otherwise "VDGBNP".
        external_rotamers : List[str]
            List of paths to external rotamer files. Default = None.
        external_templates : List[str]
            List of paths to external template files. Default = None.
        as_datalocal : bool
            Save output files to DataLocal folder. Default = True.
        pele_dir : str
            Path to PELE directory, e.g. LIG_Pele.
        """
        self.pdb = pdb_file
        self.forcefield = self._retrieve_forcefield(forcefield)
        self.charge_parametrization_method = self._check_charge_parametrization_method(
            charge_parametrization_method
        )
        self.gridres = gridres
        self.solvent = self._check_solvent(solvent, forcefield)
        self.external_templates = external_templates
        self.external_rotamers = external_rotamers
        self.as_datalocal = as_datalocal
        self.ligands = self.extract_ligands(pdb=self.pdb, gridres=self.gridres)
        self.rotamers_to_skip, self.templates_to_skip = self._check_external_files()
        self.pele_dir = pele_dir

    @classmethod
    def from_parameters(cls, parameters):
        """
        Initializes LigandParametrization from simulation parameters.

        Parameters
        -----------
        parameters : ParametersBuilder object
            Simulation parameters passed from Adaptive.simulation.

        Returns
        --------
        obj : LigandParametrization object
            LigandParametrization object initialized from simulation parameters.
        """
        obj = LigandParametrization(
            pdb_file=parameters.system,
            forcefield=parameters.forcefield,
            charge_parametrization_method=parameters.charge_parametrization_method,
            gridres=parameters.gridres,
            solvent=parameters.solvent,
            external_templates=parameters.external_template,
            external_rotamers=parameters.external_rotamers,
            as_datalocal=False,
            pele_dir=parameters.pele_dir,
        )
        return obj

    def generate(self) -> None:
        """
        Generates forcefield templates and rotamer files for ligands, then copies the one provided by the user.
        """
        for ligand in self.ligands:

            output_handler = OutputPathHandler(ligand, self.forcefield, as_datalocal=True, output_path=self.pele_dir)
            rotamer_library_path = output_handler.get_rotamer_library_path()
            impact_template_path = output_handler.get_impact_template_path()
            # solvent_template_path = output_handler.get_solvent_template_path()

            if ligand not in self.rotamers_to_skip:
                rotamer_library = RotamerLibrary(ligand)
                rotamer_library.to_file(rotamer_library_path)

            if ligand not in self.templates_to_skip:
                parameters = self.forcefield.parameterize(
                    ligand, charge_method=self.charge_parametrization_method
                )
                topology = Topology(ligand, parameters)
                impact = Impact(topology)
                impact.to_file(impact_template_path)

        self.copy_external_parameters()

    def copy_external_parameters(self) -> None:
        """
        TODO: Copy external ligand templates and rotamers specified by the user. Raised raise TemplateFileNotFound if file
        not found.
        """
        pass

    def _retrieve_forcefield(self, forcefield):
        """
        Maps forcefield YAML argument with peleffy classes.
        """
        forcefield = self.forcefields.get(forcefield.lower(), ff.OPLS2005ForceField)
        return forcefield

    def _check_charge_parametrization_method(self, method):
        """
        Checks if charge parametrization method selected by the user is supported.
        """
        if method.lower() not in constants.charge_parametrization_methods:
            raise ValueError(
                f"Invalid charge parametrization method, choose one of: {self.charge_parametrization_methods}."
            )
        else:
            return method.lower()

    def _check_solvent(self, solvent, forcefield):
        """
        Checks if solvent is compatible with the forcefield. OpenFF forcefield supports OBC solvent only.
        """
        if forcefield.lower() != "opls2005" and solvent.lower() == "vdgbnp":
            raise ValueError(
                "OpenFF support OBC solvent only. Change forcefield to 'OPLS2005' or solvent to 'OBC'."
            )
        else:
            return solvent.lower()

    @staticmethod
    def extract_ligands(pdb, gridres):
        """
        Extracts all hetero molecules in PDB and returns them as peleffy.topology.Molecule objects.
        """
        reader = PDB(pdb)
        molecules = reader.get_hetero_molecules(
            rotamer_resolution=gridres, allow_undefined_stereo=True
        )
        return molecules

    def _check_external_files(self):
        """
        Checks if any of the hetero molecules extracted from PDB has a user-defined rotamers or template file, so they
        can be skipped.
        """
        ligands = [ligand.tag for ligand in self.ligands]

        external_rotamer_residues = [
            file.split(".")[0] for file in self.external_rotamers
        ]
        rotamers_to_skip = [
            residue for residue in external_rotamer_residues if residue in ligands
        ]

        external_template_residues = [
            file.replace("z", "").upper() for file in self.external_templates
        ]
        templates_to_skip = [
            residue for residue in external_template_residues if residue in ligands
        ]

        return rotamers_to_skip, templates_to_skip
