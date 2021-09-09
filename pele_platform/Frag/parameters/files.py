"""
This module handles the assignment of file parameters for FragPELE.
"""
from pele_platform.context import context


class FragInputFiles(object):
    """
    Base class to assign the file parameters for FragPELE.
    """

    def __init__(self):
        """
        Given a Parameters object, it initializes the file parameters for
        FragPELE.

        .. todo ::
           * We need to unify all classes that prepare the parameters
             for Frag. We need an abstract class to ensure that all of them
             modify correctly the Parameters object
        """
        self.input_parameters()
        self.frag_from_core_parameters()
        self.frag_from_sdf_parameters()

    @staticmethod
    def input_parameters():
        """
        Given a Parameters object, it assigns the general parameters
        belonging to the ligand-protein input.
        """
        context.parameters.input = context.yaml_parser.frag_input
        context.parameters.frag_core_atom = context.yaml_parser.frag_core_atom

    @staticmethod
    def frag_from_core_parameters():
        """
        Given a Parameters object, it assigns the parameters that are required
        to run a FragPELE simulation with a series file.
        """
        import os
        import PPP.main as ppp

        context.parameters.core = context.yaml_parser.frag_core

        if context.yaml_parser.skip_prep is not None:
            context.parameters.skip_prep = context.yaml_parser.skip_prep
        else:
            context.parameters.skip_prep = context.parameters.simulation_params.get("skip_prep", False)

        if not context.parameters.skip_prep:
            context.parameters.core_process = \
                os.path.basename(
                    ppp.main(context.parameters.core, ".", output_pdb=["", ],
                             charge_terminals=context.yaml_parser.charge_ter,
                             no_gaps_ter=context.yaml_parser.gaps_ter)[0])
        else:
            print("Skipping preprocessing")
            context.parameters.core_process = context.parameters.core

        context.parameters.core_format = context.parameters.core.rsplit(".")[-1]

        from pathlib import Path
        import shutil

        new_file = os.path.join(os.getcwd(), Path(context.parameters.core_process).stem + "_frag.pdb")
        shutil.copy(context.parameters.core, new_file)
        context.parameters.core = new_file
        context.parameters.core_process = os.path.basename(new_file)

    @staticmethod
    def frag_from_sdf_parameters():
        """
        Given a Parameters object, it assigns the parameters that are required
        to run a FragPELE simulation with a sd file.
        """
        context.parameters.ligands = context.yaml_parser.frag_ligands
        context.parameters.frag_library = context.yaml_parser.frag_library
        context.parameters.fragment_atom = context.yaml_parser.fragment_atom
        context.parameters.frag_restart_libraries = context.yaml_parser.frag_restart_libraries
