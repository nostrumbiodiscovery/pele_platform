import os
import PlopRotTemp as plop
import pele_platform.constants.constants as cs
import pele_platform.Utilities.Helpers.helpers as hp
import pele_platform.Errors.custom_errors as ce

try:
    import subprocess32 as subprocess
except ImportError:
    import subprocess
except SyntaxError:
    import subprocess

def parametrize_miss_residues(env, resname=None):
    resname = env.residue if not resname else resname
    SPYTHON = os.path.join(cs.SCHRODINGER, "utilities/python")
    if not os.path.exists(SPYTHON):
        SPYTHON = os.path.join(cs.SCHRODINGER, "run")
    file_path = os.path.join(os.path.dirname(plop.__file__), "main.py")
    options = retrieve_options(env)
    templatedir = os.path.join(env.pele_dir, "DataLocal/Templates/OPLS2005/HeteroAtoms")
    rotamerdir = os.path.join(env.pele_dir, "DataLocal/LigandRotamerLibs")  
    mae_cahrges = True if env.mae_lig else False
    my_env = os.environ.copy()
    my_env["SCHRODINGER_PYTHONPATH"]=os.path.join(cs.SCHRODINGER, "internal/lib/python2.7/site-packages/")
    my_env["SCHRODINGER"]=cs.SCHRODINGER
    print("Running Plop")
    print("{} {} {} {} --outputname {} --templatedir {} --rotamerdir {}".format(SPYTHON, file_path, options, env.lig, resname, templatedir, rotamerdir))
    try:
        subprocess.check_output("{} {} {} {} --outputname {} --templatedir {} --rotamerdir {}".format(SPYTHON, file_path, options, env.lig, resname, templatedir, rotamerdir).split(), env=my_env)
    except subprocess.CalledProcessError: 
        raise ce.LigandPreparationError("Ligand preparation failed.\
 Plese check there are no space in the ligand atom name and that \
the inputted ligand has a valid structure")
    #hp.silentremove([syst.lig])


def retrieve_options(env):
    """
    Retrieve PlopRotTemp options from input arguments
    """

    options = []
    if env.core != -1:
        options.extend(["--core {}".format(env.core)])
    if env.mtor != 4:
        options.extend(["--mtor {}".format(env.mtor)])
    if env.n != 1000:
        options.extend(["--n {}".format(env.n)])
    if env.forcefield != "OPLS2005":
        options.extend(["--force {}".format(env.forcefield)])
    if env.mae_lig:
        options.extend(["--mae_charges"])
    if env.gridres != 10:
        options.extend(["--gridres {}".format(env.gridres)])
    return " ".join(options)

