import pele_platform.constants.constants as cs
import os
import subprocess
from time import sleep


def run_prepwizard(input):
    prepwizard = os.path.join(cs.SCHRODINGER, "utilities/prepwizard")
    output = os.path.basename(input.split(".")[0]) + "prepwizard." + input.split(".")[1]
    command = "{} {} {}".format(prepwizard, input, output)
    print(command)
    subprocess.call(command.split())
    while not os.path.exists(output):
        sleep(30)
    return output

