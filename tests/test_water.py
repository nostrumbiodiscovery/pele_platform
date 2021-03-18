from pele_platform import main
from pele_platform.constants import constants as cs
from pele_platform.Utilities.Helpers import water as wt
from pele_platform.Errors import custom_errors as ce
import shutil
import glob
import os
from . import test_adaptive as tk


test_path = os.path.join(cs.DIR, "Examples")
yaml = os.path.join(test_path, "water/input_induced.yaml")
yaml_previous_water = os.path.join(test_path, "water/input_previous_water.yaml")


def test_water(yaml=yaml):

    water_lines = [
        "HETATM 1699  OW  HOH A 202      51.000  92.000  14.000  1.00  0.00           O",
        "HETATM 1700 1HW  HOH A 202      51.757  92.586  14.000  1.00  0.00           H",
        "HETATM 1701 2HW  HOH A 202      50.243  92.586  14.000  1.00  0.00           H",
        "HETATM 1702  OW  HOH A 203      71.000  60.000  20.000  1.00  0.00           O",
        "HETATM 1703 1HW  HOH A 203      71.757  60.586  20.000  1.00  0.00           H",
        "HETATM 1704 2HW  HOH A 203      70.243  60.586  20.000  1.00  0.00           H",
        "HETATM 1705  OW  HOH A 204      82.000  86.000  74.000  1.00  0.00           O",
        "HETATM 1706 1HW  HOH A 204      82.757  86.586  74.000  1.00  0.00           H",
        "HETATM 1707 2HW  HOH A 204      81.243  86.586  74.000  1.00  0.00           H"]
 
    water_output = []

    # Function to test
    job = main.run_platform(yaml)

    # checkpoints
    output = glob.glob(os.path.join(job.pele_dir, "results/top_poses/0.1.0*"))[0]

    with open(output, "r") as file:
        lines = file.readlines()

        for line in lines:
            if line[17:21].strip() == "HOH":
                water_output.append(line.strip())

    # test
    assert water_lines == water_output

WATER_PREVIOUS_PELE = [
    '"waterSites": [{"watersToPerturb": {"links": {"ids": ["C:2091"] }}, "Box": {"radius": 6, "fixedCenter": [0.0000000000, 0.0000000000, 0.0000000000], "type": "sphericalBox"}}]',
    'WaterPerturbation::parameters',
    '"waterPerturbationFrequency": 1,'
]
WATER_PREVIOUS_ADAPTIVE = [
'"type" : "null",'
]
def test_water_with_previous_water(yaml=yaml_previous_water):
    # Function to test
    errors = []
    job = main.run_platform(yaml)
    errors = tk.check_file(job.pele_dir, "pele.conf", WATER_PREVIOUS_PELE, errors)
    errors = tk.check_file(job.pele_dir, "adaptive.conf", WATER_PREVIOUS_ADAPTIVE, errors)
    assert not errors
    
pdb = os.path.join(test_path, "gpcr/complex.pdb")
n_waters = 2
API_WATERS = [
'HETATM 4695  OW  HOH A 402      51.000  92.000  14.000  1.00  0.00           O',
'HETATM 4698  OW  HOH A 403      71.000  60.000  20.000  1.00  0.00           O'
]
def test_error_water(pdb=pdb, n_waters=n_waters):
    errors = []
    shutil.copy(pdb, ".")
    input_file = os.path.basename(pdb)
    water_object = wt.WaterIncluder([input_file], n_waters, test=True)
    try:
        water_object.run()
    except ce.NotCenterOfWaterBox:
        assert True
        return
    assert False


def test_include_water_no_ligand_API(pdb=pdb, n_waters=n_waters, water_center=[1,2,3], water_radius=18,
    water_trials=1, water_temp=1, water_constr=1, water_overlap=1, water_freq=9):
    errors = []
    shutil.copy(pdb, ".")
    input_file = os.path.basename(pdb)
    water_object = wt.WaterIncluder([input_file], n_waters, test=True,
        water_center=water_center, water_radius=water_radius, 
        water_trials=water_trials, water_temp=water_temp, water_constr=water_constr,
        water_overlap=water_overlap, water_freq=water_freq)
    water_object.run()
    errors = tk.check_file(".", input_file, API_WATERS, errors)
    #os.remove(input_file)
    assert water_object.water_line == '\n         "WaterPerturbation":\n         {\n             "watersToPerturb": { "links": { "ids": [ "A:402", "A:403" ] } },\n             "parameters":\n             {\n                 \n                 "temperature": 1,\n                 "numberOfStericTrials": 1,\n                 "overlapFactor": 1,\n                 "COMConstraintConstant": 1\n             },\n             "waterSites": [{"watersToPerturb": {"links": {"ids": ["A:402", "A:403"] }}, "Box": {"radius": 18, "fixedCenter": [1.0000000000, 2.0000000000, 3.0000000000], "type": "sphericalBox"}}]\n         }, \n'
    assert not errors
    
    
API_WATERS2 = [
'HETATM 4761  OW  HOH A 402      51.000  92.000  14.000  1.00  0.00           O',
'HETATM 4764  OW  HOH A 403      71.000  60.000  20.000  1.00  0.00           O'
]
def test_include_water_ligand_API(pdb=pdb, n_waters=n_waters, water_center=False, water_radius=False,
    water_trials=1, water_temp=1, water_constr=1, water_overlap=1, water_freq=9, ligand_residue="LIG"):
    errors = []
    shutil.copy(pdb, ".")
    input_file = os.path.basename(pdb)
    water_object = wt.WaterIncluder([input_file], n_waters, test=True,
        water_center=water_center, water_radius=water_radius, 
        water_trials=water_trials, water_temp=water_temp, water_constr=water_constr,
        water_overlap=water_overlap, water_freq=water_freq, ligand_residue=ligand_residue)
    water_object.run()
    errors = tk.check_file(".", input_file, API_WATERS2, errors)
    #os.remove(input_file)
    assert water_object.water_line == '\n         "WaterPerturbation":\n         {\n             "watersToPerturb": { "links": { "ids": [ "A:402", "A:403" ] } },\n             "parameters":\n             {\n                 \n                 "temperature": 1,\n                 "numberOfStericTrials": 1,\n                 "overlapFactor": 1,\n                 "COMConstraintConstant": 1\n             },\n             "waterSites": [{"watersToPerturb": {"links": {"ids": ["A:402", "A:403"] }}, "Box": {"radius": 6, "fixedCenter": [-87.5760388319, -7.1355193027, -64.8428620317], "type": "sphericalBox"}}]\n         }, \n'
    assert not errors
