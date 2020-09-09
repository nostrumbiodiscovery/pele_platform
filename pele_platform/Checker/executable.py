from dataclasses import dataclass
import os
import pele_platform.Errors.custom_errors as ce


@dataclass
class Executable():


    executable: str

    def is_in_path(self) -> bool:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, self.executable)
            if self._is_exe(exe_file):
                return True
        raise ce.ExecutableNotInPath(f"Executable {self.executable} not in PATH. Please do: export PATH=/path/to/folder/with/executable:$PATH")

    def _is_exe(self, fpath) -> str:
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
