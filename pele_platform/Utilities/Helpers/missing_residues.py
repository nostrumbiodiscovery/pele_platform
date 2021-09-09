import os
import pele_platform.Utilities.Helpers.plop_launcher as plop
import pele_platform.Utilities.Helpers.system_prep as sp

from pele_platform.context import context


def create_template(residue):
    output_pdb = os.path.join(context.parameters.pele_dir, "miss_residue.pdb")
    syst = sp.SystemBuilder.build_system(context.parameters.system_fix, None, residue, context.parameters.pele_dir, output=output_pdb,
                                         inputs_dir=context.parameters.inputs_dir)
    plop.parametrize_miss_residues(residue, syst.lig)
