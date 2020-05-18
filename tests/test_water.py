from pele_platform import main
from pele_platform.constants import constants as cs
import glob
import os


test_path = os.path.join(cs.DIR, "Examples")
yaml = os.path.join(test_path, "water/input_induced.yaml")


def test_water(yaml=yaml):

    water_lines = [
            "HETATM 1699  OW  HOH A 202     -36.218 -16.449  14.678  1.00  0.00           O",
            "HETATM 1700 1HW  HOH A 202     -35.620 -16.201  13.912  1.00  0.00           H",
            "HETATM 1701 2HW  HOH A 202     -36.015 -15.853  15.458  1.00  0.00           H",
            "HETATM 1702  OW  HOH A 203     -40.302 -20.787   7.671  1.00  0.00           O",
            "HETATM 1703 1HW  HOH A 203     -40.060 -19.932   7.207  1.00  0.00           H",
            "HETATM 1704 2HW  HOH A 203     -39.467 -21.320   7.826  1.00  0.00           H",
            "HETATM 1705  OW  HOH A 204     -34.970 -14.843   7.199  1.00  0.00           O",
            "HETATM 1706 1HW  HOH A 204     -35.573 -14.308   6.603  1.00  0.00           H",
            "HETATM 1707 2HW  HOH A 204     -35.522 -15.513   7.701  1.00  0.00           H"]
    
    water_output = []

    # Function to test
    job, _ = main.run_platform(yaml)

    # checkpoints
    output = glob.glob(os.path.join(job.pele_dir, "results/BestStructs/epoch0_trajectepoch0_trajectory_1.1*"))

    with open(output, "r") as file:
        lines = file.readlines()

        for line in lines:
            if line[17:21].strip() == "HOH":
                water_output.append(line.strip())

    # test
    assert water_lines == water_output
