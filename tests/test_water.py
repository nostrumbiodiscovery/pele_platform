import tests.utils
from pele_platform import main
from pele_platform.constants import constants as cs
from pele_platform.Utilities.Helpers import water as wt
from pele_platform.Errors import custom_errors as ce
import shutil
import glob
import os
import pytest


test_path = os.path.join(cs.DIR, "Examples")

WATERLIG_ARGS = os.path.join(test_path, "water/input_lig.yaml")
ALL_WATER_ARGS = os.path.join(test_path, "water/input_all.yaml")
NWATER_ARGS = os.path.join(test_path, "water/input_nwaters.yaml")

ALL_WATER_VALUES = ["WaterPerturbation::parameters", '"M:1"', '"M:2"']
WATER_VALUES = [
    "WaterPerturbation::parameters",
    '"M:1"',
]

PDB = os.path.join(test_path, "gpcr/complex.pdb")
N_WATERS = 2
API_WATERS = [
    "HETATM 4695  OW  HOH A 402      51.000  92.000  14.000  1.00  0.00           O",
    "HETATM 4698  OW  HOH A 403      71.000  60.000  20.000  1.00  0.00           O",
]


WATER_PREVIOUS_PELE = [
    '"waterSites": [{"watersToPerturb": {"links": {"ids": ["C:2091"] }}, "Box": {"radius": 6.0, "fixedCenter": [0.0000000000, 0.0000000000, 0.0000000000], "type": "sphericalBox"}}]',
    "WaterPerturbation::parameters",
    '"waterPerturbationFrequency": 1,',
]
WATER_PREVIOUS_ADAPTIVE = ['"type" : "null",']


def test_water():

    water_lines = [
        "HETATM 1699  OW  HOH A 202      51.000  92.000  14.000  1.00  0.00           O",
        "HETATM 1700 1HW  HOH A 202      51.757  92.586  14.000  1.00  0.00           H",
        "HETATM 1701 2HW  HOH A 202      50.243  92.586  14.000  1.00  0.00           H",
        "HETATM 1702  OW  HOH A 203      71.000  60.000  20.000  1.00  0.00           O",
        "HETATM 1703 1HW  HOH A 203      71.757  60.586  20.000  1.00  0.00           H",
        "HETATM 1704 2HW  HOH A 203      70.243  60.586  20.000  1.00  0.00           H",
        "HETATM 1705  OW  HOH A 204      82.000  86.000  74.000  1.00  0.00           O",
        "HETATM 1706 1HW  HOH A 204      82.757  86.586  74.000  1.00  0.00           H",
        "HETATM 1707 2HW  HOH A 204      81.243  86.586  74.000  1.00  0.00           H",
    ]

    water_output = []

    # Function to test
    yaml = os.path.join(test_path, "water/input_induced.yaml")
    job, _, _ = main.run_platform_from_yaml(yaml)

    # checkpoints
    output = glob.glob(os.path.join(job.pele_dir, "results/top_poses/0.1.0*"))[0]

    with open(output, "r") as file:
        lines = file.readlines()

        for line in lines:
            if line[17:21].strip() == "HOH":
                water_output.append(line.strip())

    # test
    assert water_lines == water_output


def test_water_with_previous_water():
    errors = []
    yaml = os.path.join(test_path, "water/input_previous_water.yaml")
    job, = main.run_platform_from_yaml(yaml)
    errors = tests.utils.check_file(job.pele_dir, "pele.conf", WATER_PREVIOUS_PELE, errors)
    errors = tests.utils.check_file(
        job.pele_dir, "adaptive.conf", WATER_PREVIOUS_ADAPTIVE, errors
    )
    assert not errors


def test_error_water(pdb=PDB, n_waters=N_WATERS):
    shutil.copy(pdb, ".")
    input_file = os.path.basename(pdb)
    water_object = wt.WaterIncluder([input_file], n_waters, test=True)
    with pytest.raises(ce.NotCenterOfWaterBox):
        water_object.run()


def test_include_water_no_ligand_API(
    pdb=PDB,
    n_waters=N_WATERS,
    water_center=[1, 2, 3],
    water_radius=18,
    water_trials=1,
    water_temp=1,
    water_constr=1,
    water_overlap=1,
    water_freq=9,
):
    errors = []
    shutil.copy(pdb, ".")
    input_file = os.path.basename(pdb)
    water_object = wt.WaterIncluder(
        [input_file],
        n_waters,
        test=True,
        water_center=water_center,
        water_radius=water_radius,
        water_trials=water_trials,
        water_temp=water_temp,
        water_constr=water_constr,
        water_overlap=water_overlap,
        water_freq=water_freq,
    )
    water_object.run()
    errors = tests.utils.check_file(".", input_file, API_WATERS, errors)
    # os.remove(input_file)
    assert (
        water_object.water_line
        == '\n         "WaterPerturbation":\n         {\n             "watersToPerturb": { "links": { "ids": [ "A:402", "A:403" ] } },\n             "parameters":\n             {\n                 \n                 "temperature": 1,\n                 "numberOfStericTrials": 1,\n                 "overlapFactor": 1,\n                 "COMConstraintConstant": 1\n             },\n             "waterSites": [{"watersToPerturb": {"links": {"ids": ["A:402", "A:403"] }}, "Box": {"radius": 18, "fixedCenter": [1.0000000000, 2.0000000000, 3.0000000000], "type": "sphericalBox"}}]\n         }, \n'
    )
    assert not errors


API_WATERS2 = [
    "HETATM 4761  OW  HOH A 402      51.000  92.000  14.000  1.00  0.00           O",
    "HETATM 4764  OW  HOH A 403      71.000  60.000  20.000  1.00  0.00           O",
]


def test_include_water_ligand_API(
    pdb=PDB,
    n_waters=N_WATERS,
    water_center=False,
    water_radius=False,
    water_trials=1,
    water_temp=1,
    water_constr=1,
    water_overlap=1,
    water_freq=9,
    ligand_residue="LIG",
):
    errors = []
    shutil.copy(pdb, ".")
    input_file = os.path.basename(pdb)
    water_object = wt.WaterIncluder(
        [input_file],
        n_waters,
        test=True,
        water_center=water_center,
        water_radius=water_radius,
        water_trials=water_trials,
        water_temp=water_temp,
        water_constr=water_constr,
        water_overlap=water_overlap,
        water_freq=water_freq,
        ligand_residue=ligand_residue,
    )
    water_object.run()
    errors = tests.utils.check_file(".", input_file, API_WATERS2, errors)

    assert (
        water_object.water_line
        == '\n         "WaterPerturbation":\n         {\n             "watersToPerturb": { "links": { "ids": [ "A:402", "A:403" ] } },\n             "parameters":\n             {\n                 \n                 "temperature": 1,\n                 "numberOfStericTrials": 1,\n                 "overlapFactor": 1,\n                 "COMConstraintConstant": 1\n             },\n             "waterSites": [{"watersToPerturb": {"links": {"ids": ["A:402", "A:403"] }}, "Box": {"radius": 6, "fixedCenter": [-87.5760388319, -7.1355193027, -64.8428620317], "type": "sphericalBox"}}]\n         }, \n'
    )
    assert not errors


def test_n_water(ext_args=NWATER_ARGS):
    job, sel, job2 = main.run_platform_from_yaml(ext_args)
    results = glob.glob(os.path.join(job.pele_dir, "results/BestStructs/*.pdb"))
    error = False
    # Result has waters
    for result in results:
        with open(result, "r") as f:
            if "HOH" not in "".join(f.readlines()):
                error = True
    # Input has no water
    with open("../pele_platform/Examples/Msm/PR_1A28_xray_-_minimized.pdb", "r") as f:
        if "HOH" in "".join(f.readlines()):
            error = True
    assert not error


def test_all_waters(ext_args=ALL_WATER_ARGS):
    errors = []
    job, = main.run_platform_from_yaml(ext_args)
    folder = job.pele_dir
    errors = tests.utils.check_file(folder, "pele.conf", ALL_WATER_VALUES, errors)
    assert not errors


def test_water_lig(ext_args=WATERLIG_ARGS):
    errors = []
    job, = main.run_platform_from_yaml(ext_args)
    folder = job.pele_dir
    errors = tests.utils.check_file(folder, "pele.conf", WATER_VALUES, errors)
    assert not errors


def test_retrieve_indices_to_track():
    """
    Run the function to check, if it correctly parses retrieved waters to a format used by Analysis.
    """
    inp = ["A:202", "B:1000", "C:1"]
    expected_out = [("A", 202), ("B", 1000), ("C", 1)]

    output = wt.WaterIncluder.retrieve_indices_to_track(inp)
    assert output == expected_out


@pytest.mark.parametrize(
    ("configuration_file", "expected_list"),
    [
        ("water_pele.conf", [("A", 202), ("A", 203), ("A", 204)]),
        ("water2_pele.conf", [("W", 1), ("A", 2405), ("A", 2406)]),
    ],
)
def test_water_ids_from_conf(configuration_file, expected_list):
    """
    Tests if we can correctly extract water indices from pele.conf.

    Parameters
    -----------
    configuration_file : str
        Path to configuration file (pele.conf).
    expected_list : list[tuple]
        List of tuples containing water molecule IDs.
    """
    from pele_platform.Utilities.Helpers.water import water_ids_from_conf

    configuration_file = os.path.join(test_path, "water", configuration_file)
    output = water_ids_from_conf(configuration_file)

    for element in expected_list:
        assert element in output
