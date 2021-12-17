from dataclasses import dataclass
import os
import pele_platform.Errors.custom_errors as ce


@dataclass
class EnvVariable:
    name: str
    variable: str
    default: str
    flag: str
    env_var: str

    def is_valid(self) -> bool:
        # Check if environment variable is correctly set
        if self.variable:
            if os.path.exists(self.variable):
                return True
        elif self.default:
            if os.path.exists(self.default):
                return True
        raise ce.EnvVariableNotFound(
            f"""{self.name} not found.\n
            1) If you have the standard installation, export the environment variable by doing:\n
            \t{self.env_var}.\n
            2) If the above does not work, define the location of {self.name} with the following flag in the input.yaml:\n
            \t{self.flag}"""
        )
