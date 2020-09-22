import matplotlib
matplotlib.use("Agg") #Avoid backend issues
import pele_platform.Checker.python_version as pv
pv.check_python_version() #Avoid python2
from argparse import ArgumentParser
import pele_platform.Utilities.Helpers.yaml_parser as yp
import pele_platform.Utilities.Helpers.launcher as lc
import pele_platform.Utilities.Parameters.pele_env as pv
import pele_platform.Checker.valid_flags as vf
import pele_platform.Errors.custom_errors as ce



def parseargs(args=[]):
    '''
    Command line parser
    '''
    parser = ArgumentParser(description='Run PELE Platform')
    parser.add_argument('input_file', type=str, help='Yaml input file')
    args = parser.parse_args(args) if args else parser.parse_args()
    return args.input_file
    
def run_platform(input_yaml):
    '''
    High level function to run a PELE job. It will:
    1) Parse the input.yaml
    2) Launch job
    3) Return job parametrs
    '''
    yaml_obj = yp.YamlParser(input_yaml, vf.VALID_FLAGS_PLATFORM)
    try:
        yaml_obj.read()
    except AttributeError:
        raise ce.WrongYamlFile("Input file: {} does not look like a correct yml file".format(input_yaml))
    job_params = lc.Launcher(yaml_obj).launch()
    return job_params

if __name__ == "__main__":
    input_file = parseargs()
    job = run_platform(input_file)
