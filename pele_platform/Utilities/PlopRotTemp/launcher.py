import os
import pele_platform.constants as cs
import pele_platform.Utilities.Helpers.helpers as hp

try:
    import subprocess32 as subprocess
except ImportError:
    import subprocess
except SyntaxError:
    import subprocess

def parametrize_miss_residues(args, env, syst):
    SPYTHON = os.path.join(cs.SCHRODINGER, "utilities/python")
    file_path = os.path.abspath(os.path.join(cs.DIR, "Utilities/PlopRotTemp/main.py"))
    options = retrieve_options(args, env)
    if args.mae_lig:
        mae_charges = True
        print("Running Plop from mae")
        print("{} {} {} {} {} {}".format(SPYTHON, file_path, options, env.mae_lig, args.residue, env.pele_dir))
        subprocess.call("{} {} {} {} {} {}".format(SPYTHON, file_path, options, env.mae_lig, args.residue, env.pele_dir).split())
        hp.silentremove([syst.system])
    else:
        mae_charges = False
        print("Running Plop from pdb")
        print("{} {} {} {} {} {}".format(SPYTHON, file_path, options, syst.lig, args.residue, env.pele_dir))
        subprocess.call("{} {} {} {} {} {}".format(SPYTHON, file_path, options, syst.lig, args.residue, env.pele_dir).split())
        hp.silentremove([syst.lig])


def retrieve_options(args, env):
    """
    Retrieve PlopRotTemp options from input arguments
    """

    options = []
    if args.core != -1:
        options.extend(["--core {}".format(args.core)])
    if args.mtor != 4:
        options.extend(["--mtor {}".format(args.mtor)])
    if args.n != 1000:
        options.extend(["--n {}".format(args.n)])
    if args.forcefield != "OPLS2005":
        options.extend(["--force {}".format(args.forcefield)])
    if args.mae_lig:
        options.extend(["--mae_charges"])
    if args.gridres != 10:
        options.extend(["--gridres {}".format(args.gridres)])
    return " ".join(options)

