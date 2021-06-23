from dataclasses import dataclass
import itertools
from packaging import version

import pele_platform.constants.constants as cs
from pele_platform.Checker import executable as ex
from pele_platform.Checker import environment_variables as en
import pele_platform.Utilities.Helpers.yaml_parser as yp
from pele_platform.Utilities.Helpers import helpers
from pele_platform.Errors.custom_errors import IncompatiblePELEVersion


EXECUTABLES_TO_CHECK = [
    "mpirun",
]


@dataclass
class Checker:
    def check_variables(self, args: yp.YamlParser) -> None:
        for env_variable in self._generate_env_variables(args):
            env_variable.is_valid()
        for executable in self._generate_executables():
            executable.is_in_path()
        # self._check_pele_version(args) TODO: Figure out where to put it (and compatibility with singularity).

    def _generate_env_variables(self, args: yp.YamlParser) -> list:
        self.env_variables = [
            en.EnvVariable(
                "pele_data",
                args.pele_data,
                cs.DEFAULT_PELE_DATA,
                "--pele_data /path/to/data/folder/",
                "export PELE_DATA=/path/to/PELE-Data/folder",
            ),
            en.EnvVariable(
                "pele_documents",
                args.pele_documents,
                cs.DEFAULT_PELE_DOCUMENTS,
                "--pele_documents /path/to/documents/folder",
                "export PELE_DOCUMENTS=/path/to/PELE-Documents/folder",
            ),
            en.EnvVariable(
                "pele_exec",
                args.pele_exec,
                cs.DEFAULT_PELE_EXEC,
                "--pele_exec /path/to/PELE_exec",
                "export PELE_EXEC=/path/to/PELE-1.X/binary",
            ),
            en.EnvVariable(
                "pele_license",
                args.pele_license,
                cs.DEFAULT_PELE_LICENSE,
                "--pele_license /path/to/licenses",
                "export PELE_LICENSE=/path/to/PELE-License/folder",
            ),
            en.EnvVariable(
                "schrodinger",
                args.schrodinger,
                cs.SCHRODINGER,
                "--schrodinger /path/to/schrodinger-20XX/",
                "export SCHRODINGER=/path/to/schrodinger-20XX/",
            ),
        ]
        return self.env_variables

    def _generate_executables(self) -> list:
        self.executables = [
            ex.Executable(executable) for executable in EXECUTABLES_TO_CHECK
        ]
        return self.executables

    def _check_pele_version(self, args: yp.YamlParser):
        """
        Checks if selected package is compatible with installed PELE version.

        Raises
        -------
        IncompatiblePELEVersion when the package is not supported.
        """
        compatibility = {
            "1.7.1": ["minimum_steps"],
            "1.7.0": ["ligand_conformations", "covalent_residue"],
        }

        current_version = helpers.get_pele_version(args.pele_exec)

        for pele_version in sorted(compatibility, reverse=True):
            if current_version < version.parse(pele_version):
                unsupported_packages = [
                    compatibility[key]
                    for key in compatibility
                    if version.parse(key) > current_version
                ]

                unsupported_packages = list(
                    itertools.chain(*unsupported_packages)
                )  # reduce the list

                for package in unsupported_packages:
                    if getattr(args, package, None):
                        raise IncompatiblePELEVersion(
                            f"Warning: {package} not supported in version PELE {current_version}."
                        )


@dataclass
class SingularityChecker(Checker):
    """
    Class to generate a valid checker for executions with singularity containers
    """

    def check_variables(self, args: yp.YamlParser) -> None:
        """
        Check all external requirements
        """
        super().check_variables(args)
        self._check_singularity_exec(args.singularity_exec)

    def _generate_env_variables(self, args: yp.YamlParser) -> list:
        """
        Define environment variables to check

        Returns
        ----------
        env_variables : list
            List of environment variables.
        """
        self.env_variables = [
            en.EnvVariable(
                "pele_license",
                args.pele_license,
                cs.DEFAULT_PELE_LICENSE,
                "--pele_license /path/to/licenses",
                "export PELE_LICENSE=/path/to/pele_license/",
            ),
            en.EnvVariable(
                "schrodinger",
                args.schrodinger,
                cs.SCHRODINGER,
                "--schrodinger /path/to/schrodinger-20XX/",
                "export SCHRODINGER=/path/to/schrodinger-20XX/",
            ),
        ]
        return self.env_variables

    def _check_singularity_exec(self, fpath) -> None:
        """
        Check if singularity container binary exists
        """
        exec_path = ex.Executable(fpath)
        exec_path.is_in_path()


def check_executable_and_env_variables(args: yp.YamlParser):
    """
    Check all external requirements are there
    before starting the simulation.
    1) Check env variables
    2) Check executables
    """
    args.singularity_exec = (
        args.singularity_exec if args.singularity_exec else cs.SINGULARITY_EXEC
    )
    checker = SingularityChecker() if args.singularity_exec else Checker()
    checker.check_variables(args)
