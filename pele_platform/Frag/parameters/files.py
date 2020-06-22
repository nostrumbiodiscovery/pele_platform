import os
import PPP.main as ppp


# Ligand-protein input
class FragInputParameters:
    def __init__(self, args):
        self.input = args.frag_input


# If using input + series_file ligands
class FragFromCoreParameters:

    def __init__(self, args):
        self.core = args.frag_core
        self.skip_prep = args.skip_prep if args.skip_prep else self.simulation_params.get("skip_prep", False)
        if not self.skip_prep:
            self.core_process = os.path.basename(ppp.main(self.core, ".", output_pdb=["", ],
                                                          charge_terminals=args.charge_ter, no_gaps_ter=args.gaps_ter)[0])
        else:
            print("Skipping preprocessing")
            self.core_process = self.core

    @property
    def core_format(self):
        return self.core.rsplit(".")[-1]


# If using input + sdf ligands
class FragFromSDFParameters:

    def __init__(self, args):
        self.ligands = args.frag_ligands


# Base class control all input file parameters
class FragInputFiles(FragInputParameters, FragFromCoreParameters, FragFromSDFParameters):

    def __init__(self, args):
        FragInputParameters.__init__(self, args)
        FragFromCoreParameters.__init__(self, args)
        FragFromSDFParameters.__init__(self, args)
