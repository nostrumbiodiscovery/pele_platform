import os
import PlopRotTemp as plop
import pele_platform.constants.constants as cs
import pele_platform.Utilities.Helpers.helpers as hp

try:
    import subprocess32 as subprocess
except ImportError:
    import subprocess
except SyntaxError:
    import subprocess

def parametrize_miss_residues(args, env, syst, resname=None):
    resname = args.residue if not resname else resname
    SPYTHON = os.path.join(cs.SCHRODINGER, "utilities/python")
    if not os.path.exists(SPYTHON):
        SPYTHON = os.path.join(cs.SCHRODINGER, "run")
    file_path = os.path.join(os.path.dirname(plop.__file__), "main.py")
    options = retrieve_options(args, env)
    templatedir = os.path.join(env.pele_dir, "DataLocal/Templates/OPLS2005/HeteroAtoms")
    rotamerdir = os.path.join(env.pele_dir, "DataLocal/LigandRotamerLibs")  
    mae_cahrges = True if args.mae_lig else False
    my_env = os.environ.copy()
    my_env["SCHRODINGER_PYTHONPATH"]=os.path.join(cs.SCHRODINGER, "internal/lib/python2.7/site-packages/")
    my_env["SCHRODINGER"]=cs.SCHRODINGER
    print("Running Plop")
    print("{} {} {} {} --outputname {} --templatedir {} --rotamerdir {}".format(SPYTHON, file_path, options, syst.lig, resname, templatedir, rotamerdir))
    subprocess.call("{} {} {} {} --outputname {} --templatedir {} --rotamerdir {}".format(SPYTHON, file_path, options, syst.lig, resname, templatedir, rotamerdir).split(), env=my_env)
    #hp.silentremove([syst.lig])


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

