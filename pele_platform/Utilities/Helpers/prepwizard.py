import pele_platform.constants.constants as cs


def run_prepwizard(input):
    prepwizard = os.path.join(cs.SCHRODINGER, "utilities/prepwizard")
    output = input.split(".")[0] + "prepwizard." + input.split(".")[1]
    command = "{} {} {}".format(prepwizard, input, output)
    subprocess.call(command.split())
    while not os.path.exists(output):
        sleep(30)
    return output

