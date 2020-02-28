import os
import shutil
import pele_platform.constants.constants as cs
import pele_platform.main as main

test_path = os.path.join(cs.DIR, "Examples")


FRAG_ARGS = [os.path.join(test_path, "frag/input.yaml")]


def test_frag(ext_args=FRAG_ARGS):
    arguments = main.parseargs_yaml(ext_args)
    arguments = main.YamlParser(arguments.input_file)
    main.Launcher(arguments).launch()



