import os
import pytest

from pele_platform.drug_design import pipelines
from pele_platform.context import context
from .utils import test_path


@pytest.mark.parametrize(
    ("package", "yaml_file"),
    [
        ("pipelines.SiteFinder", os.path.join("site_finder", "input_folder.yaml")),
        ("pipelines.PPI", os.path.join("PPI", "input.yaml")),
        ("pipelines.GPCR", os.path.join("gpcr", "input.yaml")),
        ("pipelines.OutIn", os.path.join("out_in", "input.yaml")),
        ("pipelines.InducedFitFast", os.path.join("induced_fit", "input_fast.yaml")),
        (
            "pipelines.InducedFitExhaustive",
            os.path.join("induced_fit", "input_exhaustive.yaml"),
        ),
        ("pipelines.CovalentDocking", os.path.join("covalent_docking", "input1.yaml")),
    ],
)
def test_api_packages(package, yaml_file):
    """
    Builds parameters and runs predefined package (Pipeline) for API.
    """
    yaml_file = os.path.join(test_path, yaml_file)

    context.build_parameters(yaml_file=yaml_file)
    eval(f"{package}()").run()
    context.reset_parameters()
