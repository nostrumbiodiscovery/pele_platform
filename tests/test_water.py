from pele_platform import main
from pele_platform.constants import constants as cs
import glob
import os


test_path = os.path.join(cs.DIR, "Examples")
yaml = os.path.join(test_path, "water/input_induced.yaml")


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
    output = glob.glob(os.path.join(job.pele_dir, "results/BestStructs/epoch0_trajectory_1.0*"))[0]

    with open(output, "r") as file:
        lines = file.readlines()

        for line in lines:
            if line[17:21].strip() == "HOH":
                water_output.append(line.strip())

    # test
    assert water_lines == water_output
