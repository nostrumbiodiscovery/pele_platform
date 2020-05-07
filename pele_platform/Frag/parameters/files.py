import os
import PPP.main as ppp

#Ligand-protein input
class FragInputParameters():
    def __init__(self, args):
        self.input = args.frag_input

#If using input + serie_file ligands
class FragFromCoreParameters():

    def __init__(self, args):
        self.core = args.frag_core
        self.core_process = os.path.basename(ppp.main(self.core, ".", output_pdb=["" , ],
            charge_terminals=args.charge_ter, no_gaps_ter=args.gaps_ter)[0])
        self.core_format = args.frag_core.rsplit(".")[-1]

#If using input + sdf ligands
class FragFromSDFParameters():

    def __init__(self, args):
        self.ligands = args.frag_ligands

#Base class control all input file parameters
class FragInputFiles(FragInputParameters,
    FragFromCoreParameters, FragFromSDFParameters):

    def __init__(self, args):
        FragInputParameters.__init__(self, args)
        FragFromCoreParameters.__init__(self, args)
        FragFromSDFParameters.__init__(self, args)
