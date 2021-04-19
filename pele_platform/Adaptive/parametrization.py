import os
import shutil
import subprocess
import warnings

from peleffy import solvent
from peleffy.utils.input import PDB
from peleffy.utils import OutputPathHandler
from peleffy.topology import RotamerLibrary, Topology
from peleffy import forcefield as ff
from peleffy.template import Impact

from pele_platform.Errors import custom_errors
from pele_platform.constants import constants


class Parametrization:
    # Fallback paths, if only external templates available
    OPLS_IMPACT_TEMPLATE_PATH = "DataLocal/Templates/OPLS2005/HeteroAtoms/"
    OFF_IMPACT_TEMPLATE_PATH = "DataLocal/Templates/OpenFF/Parsley/"
    ROTAMER_LIBRARY_PATH = "DataLocal/LigandRotamerLibs/"

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

    # Mapping between args.solvent and peleffy models
    solvents = {
        "obc2": solvent.OBC2,
        "opls_obc": solvent.OPLSOBC,
        "vdgbnp": None,
    }

    # Available methods of charge parametrization
    charge_parametrization_methods = ["am1bcc", "gasteiger", "opls2005"]

    def __init__(
        self,
        pdb_file,
        forcefield="OPLS2005",
        charge_parametrization_method="am1bcc",
        gridres=10,
        solvent="OBC",
        external_templates=None,
        external_rotamers=None,
        as_datalocal=False,
        pele_dir=None,
        exclude_terminal_rotamers=True,
        ligand_core_constraints=None,
        ligand_resname=None,
    ):
        """
        Initializes Parametrization to generate template and rotamer files.

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
            Save output files to DataLocal folder. Default = False.
        pele_dir : str
            Path to PELE directory, e.g. LIG_Pele.
        exclude_terminal_rotamers : bool
            Toggle to exclude terminal rotamers.
        ligand_core_constraints : List[str]
            List of PDB atom names to be constrained as core. Default = None
        ligand_resname : str
            Residue name of the ligand. Default = None.
        """
        self.pdb = pdb_file
        self.forcefield = self._retrieve_forcefield(forcefield)
        self.charge_parametrization_method = self._check_charge_parametrization_method(
            charge_parametrization_method, forcefield
        )
        self.gridres = gridres
        self.solvent = self._retrieve_solvent_model(solvent, forcefield)
        self.external_templates = (
            external_templates if external_templates is not None else list()
        )
        self.external_rotamers = (
            external_rotamers if external_rotamers is not None else list()
        )
        self.as_datalocal = as_datalocal
        self.hetero_molecules = self.extract_ligands(
            pdb_file=self.pdb,
            gridres=self.gridres,
            exclude_terminal_rotamers=exclude_terminal_rotamers,
            ligand_resname=ligand_resname,
            ligand_core_constraints=ligand_core_constraints,
        )
        self.rotamers_to_skip, self.templates_to_skip = self._check_external_files()
        self.pele_dir = pele_dir

    @classmethod
    def from_parameters(cls, parameters):
        """
        Initializes Parametrization from simulation parameters (e.g. passed from Adaptive.simulation).

        Parameters
        -----------
        parameters : ParametersBuilder object
            Simulation parameters passed from Adaptive.simulation.

        Returns
        --------
        obj : Parametrization object
            Parametrization object initialized from simulation parameters.
        """
        obj = Parametrization(
            pdb_file=parameters.system,
            forcefield=parameters.forcefield,
            charge_parametrization_method=parameters.charge_parametrization_method,
            gridres=parameters.gridres,
            solvent=parameters.solvent,
            external_templates=parameters.external_template,
            external_rotamers=parameters.external_rotamers,
            as_datalocal=True,
            pele_dir=parameters.pele_dir,
            exclude_terminal_rotamers=parameters.exclude_terminal_rotamers,
            ligand_core_constraints=parameters.core,
            ligand_resname=parameters.resname,
        )
        return obj

    @staticmethod
    def extract_ligands(
        pdb_file,
        gridres,
        exclude_terminal_rotamers=True,
        ligand_core_constraints=None,
        ligand_resname=None,
    ):
        """
        Extracts all hetero molecules in PDB and returns them as peleffy.topology.Molecule objects.

        Parameters
        -------------
        pdb_file : str
            Path to PDB file.
        gridres : int
            Resolution of the rotamers when sampling.
        exclude_terminal_rotamers : bool
            Toggle to exclude terminal rotamers. Default = True.
        ligand_core_constraints : List[str]
            List of PDB atom names to be constrained as core. Default = None
        ligand_resname : str
            Residue name of the ligand. Default = None.

        Returns
        ---------
        unique_molecules : List[peleffy.topology.Molecule]
            List of hetero molecules extracted from the PDB file, without any duplicates, water molecules or single atom
            anions (e.g. Cl-, F-, etc.).
        """
        reader = PDB(pdb_file)
        molecules = reader.get_hetero_molecules(
            rotamer_resolution=gridres,
            allow_undefined_stereo=True,
            exclude_terminal_rotamers=exclude_terminal_rotamers,
            ligand_core_constraints=ligand_core_constraints,
            ligand_resname=ligand_resname,
        )

        # Filter out water molecules and single ions.
        to_remove = []
        for molecule in molecules:

            if molecule.tag == "HOH":
                to_remove.append(molecule)
            if (
                molecule.tag.upper().strip() in constants.ions9
                and len(molecule.get_pdb_atom_names()) == 1
            ):
                to_remove.append(molecule)

        molecules = [molecule for molecule in molecules if molecule not in to_remove]

        # Remove duplicates residue names
        unique_molecules = []

        for molecule in molecules:
            if molecule.tag not in [unique.tag for unique in unique_molecules]:
                unique_molecules.append(molecule)

        return unique_molecules

    def _copy_external_parameters(self, rotamer_path, template_paths):
        """
        Copy external ligand templates and rotamers specified by the user to both OPLS2005 and OpenFF directories, since
        we do not know which FF was used to generate them. Raised raise TemplateFileNotFound RotamersFileNotFound error,
        if file not found.

        Parameters
        -----------
        rotamer_path : str
            Path to rotamers directory in pele_dir.
        template_paths : List[str]
            Paths to template directories in pele_dir (OPLS2005 and OpenFF).
        """
        if self.external_templates:
            for file in self.external_templates:
                try:
                    # Copy into both OPLS2005 and OpenFF template directories, since we don't know what's in the file.
                    shutil.copy(file, template_paths[0])
                    shutil.copy(file, template_paths[-1])
                    print(
                        "Copied external template files:",
                        ", ".join(self.external_templates),
                    )
                except IOError:
                    raise custom_errors.TemplateFileNotFound(
                        f"Could not locate {file} file. Please double-check the path."
                    )

        if self.external_rotamers:
            for file in self.external_rotamers:
                try:
                    shutil.copy(file, rotamer_path)
                    print(
                        "Copied external rotamer files:",
                        ", ".join(self.external_rotamers),
                    )
                except IOError:
                    raise custom_errors.RotamersFileNotFound(
                        f"Could not locate {file} file. Please double-check the path."
                    )

    @staticmethod
    def _check_solvent(solvent, forcefield):
        """
        Checks if solvent is compatible with the forcefield. OpenFF forcefield supports OBC solvent only.

        Parameters
        -----------
        solvent : str
            Solvent selected by the user.
        forcefield : str
            Forcefield selected by the user.

        Returns
        --------
        solvent : str
            Solvent string, if it is compatible. Otherwise raises ValueError.
        """
        if forcefield.lower() != "opls2005" and solvent.lower() == "vdgbnp":
            raise ValueError(
                "OpenFF supports OBC solvent only. Change forcefield to 'OPLS2005' or solvent to 'OBC2' or 'OBC1'."
            )
        else:
            return solvent.lower()

    def _retrieve_forcefield(self, forcefield):
        """
        Maps forcefield YAML argument with peleffy classes.

        Parameters
        -------------
        forcefield : str
            Forcefield defined by the user in args.forcefield.

        Returns
        ---------
        forcefield : peleffy.forcefield.X
            Forcefield object from peleffy, where X depends on the forcefield selected by the user.
        """
        forcefield = self.forcefields.get(forcefield.lower(), ff.OPLS2005ForceField)
        return forcefield

    def _check_charge_parametrization_method(self, method, forcefield):
        """
        Checks if charge parametrization method selected by the user is supported and compatible with the forcefield.

        Parameters
        ------------
        method : str
            Method of charge parametrization selected by the user.
        forcefield : str
            Forcefield selected by the user.

        Returns
        --------
        method : str
            Method of charge parametrization compatible with the forcefield.
        """
        if method.lower() not in self.charge_parametrization_methods:
            raise ValueError(
                f"Invalid charge parametrization method, choose one of: {self.charge_parametrization_methods}."
            )

        # Override user-defined method, if not compatible with OPLS2005 force field.
        if forcefield.lower() == "opls2005" and method.lower() != "opls":
            print(
                f"Charge parametrization method {method} incompatible with {forcefield}. Defaulting to OPLS2005."
            )
            method = "opls2005"

        return method.lower()

    def _check_external_files(self):
        """
        Checks if any of the hetero molecules extracted from PDB has a user-defined rotamers or template file, so they
        can be skipped.

        Returns
        --------
        rotamers_to_skip : List[str]
            List of hetero molecules for which the rotamers have been supplied in an external file.
        templates_to_skip : List[str]
            List of hetero molecules for which the templates have been supplied in an external file.
        """
        ligands = [ligand.tag.strip() for ligand in self.hetero_molecules]

        if self.external_rotamers:
            # Get residue names from rotamer files.
            external_rotamer_residues = [
                os.path.basename(file).split(".")[0] for file in self.external_rotamers
            ]
            rotamers_to_skip = [
                residue for residue in external_rotamer_residues if residue in ligands
            ]
        else:
            rotamers_to_skip = list()

        if self.external_templates:
            all_templates = self.external_templates + constants.in_pele_data
        else:
            all_templates = constants.in_pele_data

        external_template_residues = [
            os.path.basename(file).rstrip("z").upper() for file in all_templates
        ]
        templates_to_skip = [
            residue for residue in external_template_residues if residue in ligands
        ]

        return rotamers_to_skip, templates_to_skip

    def generate_ligand_parameters(self):
        """
        Generates forcefield templates and rotamer files for ligands, then copies the one provided by the user.
        """
        rotamer_library_path, impact_template_paths = None, None

        for molecule in self.hetero_molecules:

            output_handler = OutputPathHandler(
                molecule,
                self.forcefield,
                as_datalocal=self.as_datalocal,
                output_path=self.pele_dir,
            )
            rotamer_library_path = output_handler.get_rotamer_library_path()
            impact_template_path = output_handler.get_impact_template_path()
            solvent_template_path = output_handler.get_solvent_template_path()

            if molecule.tag.strip() not in self.rotamers_to_skip:
                rotamer_library = RotamerLibrary(molecule)
                rotamer_library.to_file(rotamer_library_path)

            if molecule.tag.strip() not in self.templates_to_skip:
                # Try to parametrize with OPLS if OpenFF fails
                try:
                    parameters = self.forcefield.parameterize(
                        molecule, charge_method=self.charge_parametrization_method
                    )
                except (subprocess.CalledProcessError, TypeError):
                    default = "OPLS2005"
                    fallback_forcefield = self._retrieve_forcefield(default)

                    try:
                        parameters = fallback_forcefield.parameterize(molecule)
                        warnings.warn(
                            f"Could not parametrize residue {molecule.tag.strip()} with the selected forcefield. "
                            f"Parametrized with {default} instead. "
                        )
                    except subprocess.CalledProcessError as e:
                        raise custom_errors.LigandPreparationError(
                            f"Could not parametrize {molecule.tag.strip()}. The error was {e}."
                        )

                topology = Topology(molecule, parameters)
                impact = Impact(topology)
                impact.to_file(impact_template_path)

                impact_template_paths = [impact_template_path]

                print(f"Parametrized {molecule.tag.strip()}.")

                if self.solvent:
                    solvent_parameters = self.solvent(topology)
                    solvent_parameters.to_file(solvent_template_path)

        if not rotamer_library_path:
            rotamer_library_path = os.path.join(
                self.pele_dir, self.ROTAMER_LIBRARY_PATH
            )

        if not impact_template_paths:
            impact_template_paths = [
                os.path.join(self.pele_dir, self.OPLS_IMPACT_TEMPLATE_PATH),
                os.path.join(self.pele_dir, self.OFF_IMPACT_TEMPLATE_PATH),
            ]

        self._copy_external_parameters(
            os.path.dirname(rotamer_library_path),
            [os.path.dirname(path) for path in impact_template_paths],
        )

    def _retrieve_solvent_model(self, solvent, forcefield):
        """
        Checks solvent compatibility with the forcefield and returns the solvent class from peleffy.

        Parameters
        -----------
        solvent : str
            Solvent defined by the user in YAML.
        forcefield : str
            Forcefield defined by the user in YAML.

        Returns
        --------
        solvent_class : peleffy.solvent.X
            Solvent object from peleffy, X depends on the selected solvent type.
        """
        if not solvent:
            return None
        else:
            checked_solvent = self._check_solvent(solvent, forcefield)
            solvent_class = self.solvents.get(checked_solvent.lower(), None)
            return solvent_class
