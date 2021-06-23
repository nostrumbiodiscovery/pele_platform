import os
import pele_platform.Utilities.Helpers.plop_launcher as plop
import pele_platform.Utilities.Helpers.system_prep as sp


def create_template(env, residue):
    output_pdb = os.path.join(env.pele_dir, "miss_residue.pdb")
    syst = sp.SystemBuilder.build_system(env.system_fix, None, residue, env.pele_dir, output=output_pdb,
                                         inputs_dir=env.inputs_dir)
    plop.parametrize_miss_residues(env, residue, syst.lig)
