#!/usr/bin/env python
'''
Created in July, 2014
Author: Christoph Grebner
E-Mail: christoph.grebner@astrazeneca.com

program to calculate ANM and PCA for a given set of proteins in PDB format
uses the ProDy library 

'''
import prody
import glob
import argparse
import numpy as np


class PCA_ANM_Analysis:
    '''
    class definition for PCA and ANM analysis using ProDy library
    calculate the PCA of a given set of structures and the ANM of the reference structure
    '''
    def __init__(self, selection="calpha"):
        #Initialization of variables with default values
        self.selection=selection
        self.vmd=False
        self.compare=False
        self.ref="none"
        self.debug=False
        self.pdbarg="all"

        self.pdbs = []
 
    def do_Auto_analysis(self):
        pdbs = self.getPDBs()
        ensemble = self.createEnsemble(pdbs)
        pca = self.calcPCA(ensemble)
        anm = self.calcANM(ensemble.getConformation(0))
        if self.compare:
            self.compare_modes(anm,pca)
    
    def printSettings(self):
        '''
        Print all settings for debugging purposes
        '''
        print("Settings in PCA_ANM_Analysis")
        print("selection, VMD, compare, ref, pdb")
        print(self.selection, self.vmd, self.compare, self.ref, self.pdbarg)


    def set_default_based_on_argparse(self, selection, vmd, compare, ref, debug, pdb):
        '''
        set arguments from command line based on argparser
        '''
        if selection in ["calpha", "backbone", "all"]:
            self.selection=selection
        else:
            print("Invalid value for selection:", selection)
            print("Possible Values for selection are: calpha, backbone, all")
            exit()
        self.vmd = vmd
        self.compare=compare
        self.ref=ref
        self.debug=debug
        self.pdbarg=pdb
        
        #set prody level of output
        if self.debug:
            prody.confProDy(verbosity="debug")
        else:
            prody.confProDy(verbosity="info")
            
        if self.debug:
            print("***DEBUG***")
            self.printSettings()
            print("***DEBUG-END***")
    
    def getPDBs(self):
        '''
        Evaluate information f#from PCA_analysis_GUI import *
from PCA_dPCA import *
rom command line
        Get PDBs and return array
        '''
        if self.pdbarg == "all":
            pdbs = glob.glob('*.pdb')
        else:
            pdbs = self.pdbarg.split()
            
        if self.debug:
            print("***DEBUG*** The chosen PDBs are:", pdbs)
        self.pdbs = pdbs    
        return pdbs
        
        ''' 
        Possibilities
        all PDBs --> simple: pdbs = glob.glob('*.pdb')
        list with PDBs --> simple: pdbs = self.args.pdb.split()
        trajectory: check PDB if more than one structure and set flag: will be handled in createEnsemble --> done here
        DCD: ? not important
        '''

        
    def createEnsemble(self, pdbs):
        '''
        Create a prody ensemble based on getPDBs return
        Take into account, that system can be prepared or not
        and take gaps and multiple chains into account
        take single files and trajectory file into account
        
        pdbs: list-array with pdb-filenames
        
        ToDo: 
            - set reference depending on longest chain with no gaps
            - check for duplicate chains and only select this one
        '''
        print("Create Ensemble")
        ref_chids = []
        ensemble_ref_title = "Default"
        
        
        #set the reference structure, if not chosen by user, the first Frame or PDB is taken
        #remove reference Frame/PDB from list
        if self.ref == "none":
            self.ref = pdbs[0]
            pdbs.pop(0)
        else:
            if self.ref in pdbs:
                if self.debug:
                    print("***DEBUG*** Found Ref:", self.ref)
                    print("***DEBUG*** Index Ref", pdbs.index(self.ref))
                    print("Removing reference from pdblist")
                pdbs.pop(pdbs.index(self.ref))
   
        if self.debug:
            print("***DEBUG*** Reference is",self.ref)
            
        try: 
            f = open(self.ref)
        except IOError as e:
            print("I/O error({0}): {1} \"{2}\"".format(e.errno, e.strerror, self.ref))
            exit()
        
        
        
        #open reference file and check for HID, HIE; if found, replace by HIS
        try:
            pdbinfile = open(self.ref)   
        except IOError as e:
            print(self.ref)
            print("I/O error({0}): {1}".format(e.errno, e.strerror))
            exit()
        
        pdbfiledata = pdbinfile.read()
        for line in pdbinfile:
            if line[17:20] == "HID" or line[17:20] == "HIE":
                print("HID or HIE found in pdb-structure. The names will be changed to HIS in the input file")
                break
        pdbinfile.close()
        
        pdbfiledata = pdbfiledata.replace("HID","HIS")
        pdbfiledata = pdbfiledata.replace("HIE","HIS")
        
        print("Warning!!! Rewriting reference file {}".format(self.ref))
        pdbfile = open(self.ref,'w')
        pdbfile.write(pdbfiledata)
        pdbfile.close()
        
        #set reference and create first structure  
        if self.selection == "all":
            ref_structure = prody.parsePDB(self.ref)
        else:
            ref_structure = prody.parsePDB(self.ref, subset=self.selection)
    

        
        ref_selection = ref_structure.select(self.selection)
        reference_hierview = ref_structure.getHierView()
        
        print("Found", reference_hierview.numChains(), "Chain(s) in", reference_hierview._atoms.getTitle())
        ensemble_ref_title = reference_hierview._atoms.getTitle()
        print(reference_hierview[0])
        
        '''
        at the moment all chains are taken
        --> check for duplicate chains and only take on of those
        '''
        
        for chain in reference_hierview:
            ref_chids.append(chain.getChid())
        if self.debug:
            print("***DEBUG***", ref_chids)
        
        reference_chains = [reference_hierview[chid] for chid in ref_chids] 
        
        if self.debug:
            print("***DEBUG***", reference_chains)
            
        
        
        ref_chain = reference_chains[0]
        for i in range (1, len(reference_chains), 1):
            ref_chain = ref_chain + reference_chains[i]
        if self.debug:
            print("***DEBUG***", ref_chain)
        # save globally
        self.ref_chain = ref_chain
        #Create Ensemble of structures
        ensemble = prody.PDBEnsemble(ensemble_ref_title)
        ensemble.setAtoms(ref_chain)
        ensemble.setCoords(ref_chain)
        #Set ref_structure as first coordinate set
        ensemble.addCoordset(ref_structure)
        
        if self.debug:
            print("***DEBUG***", ref_chain.getResnames())
            print("***DEBUG***", ref_chain.getResnums())
        
        unmapped = []
        # map remaining structures to reference chain and add to ensemble if mapped
        for pdb in pdbs:
             if self.debug:
                 print("***DEBUG*** Processing ", pdb)
             if self.selection == "all":
                 #structure = prody.parsePDB(self.ref)
                 structure = prody.parsePDB(pdb)
             else:
                 structure = prody.parsePDB(pdb, subset=self.selection)
             
             atommaps = []
             for reference_chain in reference_chains:
                 # Map current PDB file to the reference chain
                 mappings = prody.mapOntoChain(structure, reference_chain,
                                         seqid=90,
                                         coverage=50,
                                         subset=self.selection)
                 #print(mappings, len(mappings)
                 if len(mappings) == 0:
                     print('Failed to map', pdb)
                     break
                 atommaps.append(mappings[0][0])
                 # Make sure all chains are mapped
             if len(atommaps) != len(reference_chains):
                 unmapped.append(pdb)
                 continue
             #print(atommaps[0]
             atommap = atommaps[0]
             for i in range (1, len(reference_chains), 1):
                  atommap = atommap + atommaps[i]
             # Add the atommap (mapped coordinates) to the ensemble
             # Note that some structures do not completely map (missing residues)
             # so we pass weights (1 for mapped atoms, 0 for unmapped atoms)
             ensemble.addCoordset(atommap, weights=atommap.getFlags('mapped'))
     
        if self.debug:
            print("***DEBUG*** Number of conformations:", ensemble.numConfs())
            print("***DEBUG*** ", ensemble)
            print("***DEBUG*** Unmapped structures:", unmapped)
        ensemble.iterpose()

        #print(ensemble.getTitle()
        return ensemble
    
    def setSelections(self, ensemble):
        '''
        put selections also into key arguments
        '''        
        
        custom_selection = False
        
        if custom_selection == False:
            self.selection_ref_structure = self.ref_chain
            return ensemble
        else:
            self.selection2="name CA"
            structure =  ensemble.getAtoms().select(self.selection2)
            sel_string = structure.getSelstr()
            #selec = structure.getSelstr().split()
            selec = structure.getIndices()
            if self.debug:
                print("***DEBUG***", selec)
                print("***DEBUG***", len(selec))
        
            #create new ensemble containing only the selected atoms        
            ensemble2 = prody.PDBEnsemble(ensemble.getTitle())
            ensemble2.setAtoms(structure)
            ensemble2.setCoords(structure.getCoords())
            self.selection_ref_structure = structure
            
            #ensemble2.addCoordset(structure)                
            for coords in ensemble.iterCoordsets():
                new_coords = np.zeros((ensemble2.numSelected(), 3))
                count = 0
                for coord in range(len(coords)):
                    if coord in selec:
                        new_coords[count] = coords[coord]
                        #print(count, new_coords[count]
                        count+=1
                ensemble2.addCoordset(new_coords) 
                
            ensemble2.getConformation(0).setLabel(ensemble2.getTitle())    

            if self.debug:
                print("***DEBUG*** Ensemble before selections: ", repr(ensemble))
                print("***DEBUG*** Ensemble after selections: ", repr(ensemble2))
            
            return ensemble2
    
    def setWeights(self, ensemble):
        '''
        create a custom selection (e.g. only specific residues) and set the weights for the remaining atoms
        to 0
        --> They will not contribute to the PCA
        the original ensemble is not changed!
        '''
        print("Set weights")
        self.selection2="(backbone)"
        structure =  ensemble.getAtoms().select(self.selection2)
        sel_string = structure.getSelstr()
        #selec = structure.getSelstr().split()
        selec = structure.getIndices()
        #set weights for atoms outside selection to zero!
        weights =  ensemble.getWeights()
        for i in range(1,len(weights)):
            for j in range(len(weights[i])):               
                #set atoms in selection to 1 rest to 0
                if j not in selec:
                    weights[i][j] = 0
                #print(j, weights[i][j]
        
        ensemble.setWeights(weights)

        return ensemble
        
    def calcPCA(self, ensemble):
        '''
        calcPCA:
        #ensemble: prody ensmeble with structure information
        
        calculate PCA for a set of structures    
        
        return: prody.pca object
        '''
        print("Calculate PCA")
             
        PCAname = ensemble.getTitle()
        pca = prody.PCA(PCAname)
        pca.buildCovariance(ensemble)

        print("PCA")
        pca.calcModes()
        print(repr(pca))
        
        outputname = PCAname + "_pca_modes.nmd"
                    
        prody.writeNMD(outputname, pca[:10], self.selection_ref_structure)
        
        if self.vmd == True:
            prody.viewNMDinVMD(outputname)
        print("PCA is saved in:", outputname)
        return pca, outputname
        
    def calcANM(self, structure):
        '''
        calcANM:
        #structure: prody PDB-structure
        
        calculate ANM for one specific structure    
        
        return: prody.anm object
        '''
        print("Calculate ANM")
        ANMname = structure.getLabel()
        anm = prody.ANM(ANMname)
        anm.buildHessian(structure, cutoff=15.0)
        anm.calcModes()
        # write out the three lowest modes as NMD file (visualize with NMWizard in VMD
        outputname = ANMname + "_anm_modes.nmd"
        prody.writeNMD(outputname, anm[:10], self.selection_ref_structure)
        if self.vmd == True:
            prody.viewNMDinVMD(outputname)
        print("ANM is saved in:", outputname)
        return anm

    def compare_modes(self, modes1, modes2):
        '''
        Compare the calculated ANM and PCA
        '''
        
        print(self.print_Eigen_Val_Vec(modes1))
        print(self.print_Eigen_Val_Vec(modes2))
        
        prody.printOverlapTable(modes2[:3], modes1[:3]) # Top 3 PCs vs slowest 3 ANM modes
        
        for mode in modes2[:10]:
            var = prody.calcFractVariance(mode).round(3)
            print('{0:s}  % variance = {1:.2f}'.format(mode, var))
            coll = prody.calcCollectivity(mode)
            print('{0:s}  collectivity = {1:.2f}'.format(mode, coll))
    
        for mode in modes1[:10]:    # Print ANM mode collectivity
            var = prody.calcFractVariance(mode).round(3)
            print('{0:s}  % variance = {1:.2f}'.format(mode, var))
            coll = prody.calcCollectivity(mode)
            print('{0:s}  collectivity = {1:.2f}'.format(mode, coll))
        
    def print_Eigen_Val_Vec(self,mode):
        print("Eigenvalues of", mode.getModel())
        print(mode.getEigvals().round(3))
        print("Eigenvectors of", mode.getModel())
        print(mode.getEigvecs().round(3))

    def print_Match(self, match):
       print('Chain 1     : {}'.format(match[0]))
       print('Chain 2     : {}'.format(match[1]))
       print('Length      : {}'.format(len(match[0])))
       print('Seq identity: {}'.format(match[2]))
       print('Seq overlap : {}'.format(match[3]))
       print('RMSD        : {}\n'.format(calcRMSD(match[0], match[1])))

    def print_modes(self, modes):
        print("Printing modes to file")
        import matplotlib.pyplot as plt
        residues = self.ref_chain.getResnums()
        print_offset = residues[0]
        plt.rcParams.update({'font.size': 14})

        for count in range(len(modes)):
            #array with lenghts of vectors
            a3d = (modes[count].getArrayNx3()**2).sum(axis=1)**0.5
            show = plt.plot(a3d[:])
            plt.xlabel('Residue index')
            plt.ylabel('Lenght of fluctuation vector')         
            locs,labels = plt.xticks()
            new_labels = ['%d' % (a+print_offset) for a in locs]
            plt.xticks(locs,new_labels)
            plt.title(str(modes)+" mode "+str(count+1))
            plt.savefig(str(modes)+"_mode_"+str(count)+".pdf")
            plt.clf()   
        
        for count in range(6):
            #array with lenghtes of vectors
            a3d = (modes[count].getArrayNx3()**2).sum(axis=1)**0.5
            show = plt.plot(a3d[:], label=("Mode "+str(count+1)))

        plt.xlabel('Residue index')
        plt.ylabel('Lenght of fluctuation vector')
        
        locs,labels = plt.xticks()
        new_labels = ['%d' % (a+print_offset) for a in locs]
        plt.xticks(locs,new_labels)
        lgd = plt.legend(loc='center right', bbox_to_anchor=(1.3,0.5))
        plt.savefig(str(modes)+"_mode_1_to_6.pdf", bbox_extra_artists=(lgd, ), bbox_inches='tight')
        
        print("Printing modes to file finished")

          
    def print_diff(self, pca_modes, anm_modes):
        print("here I am in DIFF")
        import matplotlib.pyplot as plt
        from matplotlib.path import Path
        import matplotlib.patches as patches
        residues = self.ref_chain.getResnums()
        print_offset = residues[0]
        #Paper01:PCA-X-Ray:
        print_offset = 738
        plt.rcParams.update({'font.size': 24}) 
        plt.rcParams['xtick.major.pad'] = 10
        plt.rcParams['ytick.major.pad'] = 10
        
        a3d_pca = ((pca_modes[0].getArrayNx3())**2).sum(axis=1)**0.5
        a3d_anm = ((anm_modes[0].getArrayNx3())**2).sum(axis=1)**0.5
        a3d_all = a3d_pca - a3d_anm
        for count in range(1,6):
            #array with lenghts of vectors
            a3d_pca = ((pca_modes[count].getArrayNx3())**2).sum(axis=1)**0.5
            a3d_anm = ((anm_modes[count].getArrayNx3())**2).sum(axis=1)**0.5
            a3d_all = a3d_all + a3d_pca - a3d_anm
            
            
        show = plt.plot(a3d_all[:], label=("Difference"), lw=1.5)
        plt.xlabel('Residue index')
        plt.ylabel('|PCA| - |ANM|')         
        locs,labels = plt.xticks()
        new_labels = ['%d' % (a+print_offset) for a in locs]
        plt.xticks(locs,new_labels)
        lgd = plt.legend(loc='center right', bbox_to_anchor=(1.65,0.5))
        ##add boxes
        codes = [Path.MOVETO,
                  Path.LINETO,
                  Path.LINETO,
                  Path.LINETO,
                  Path.CLOSEPOLY,
                  ]
        ### light green
        verts = [
                 (95, -2.0), # left, bottom
                 (95, 2.0), # left, top
                 (115, 2.0), # right, top
                 (115, -2.0), # right, bottom
                 (0., 0.), # ignored
                 ]
        path = Path(verts, codes)
        #add box around area of interest
        ax = plt.subplot()
        patch = patches.PathPatch(path, facecolor='#CCFFCC', lw=1)
        #patch = patches.PathPatch(path, facecolor='lightgreen', lw=1)
        ax.add_patch(patch)
        plt.axhline(y=0.0,xmin=0,xmax=3,c="black",linewidth=1.0,zorder=0)
        plt.savefig(str(pca_modes)+"_all_DIFF.pdf", bbox_extra_artists=(lgd, ), bbox_inches='tight')
        plt.clf()   
        
        
        for count in range(6):
            #array with lenghtes of vectors         
            a3d = ((pca_modes[count].getArrayNx3())**2).sum(axis=1)**0.5 - ((anm_modes[count].getArrayNx3())**2).sum(axis=1)**0.5         
            show = plt.plot(a3d[:], label=("Mode "+str(count+1)), lw=1.5)
        plt.xlabel('Residue index')
        plt.ylabel('|PCA| - |ANM|')   
        locs,labels = plt.xticks()
        new_labels = ['%d' % (a+print_offset) for a in locs]
        plt.xticks(locs,new_labels)
        lgd = plt.legend(loc='center right', bbox_to_anchor=(1.65,0.5))
        
        ##add boxes
        codes = [Path.MOVETO,
                  Path.LINETO,
                  Path.LINETO,
                  Path.LINETO,
                  Path.CLOSEPOLY,
                  ]
        ### light green
        verts = [
                 (95, -0.7), # left, bottom
                 (95, 0.7), # left, top
                 (115, 0.7), # right, top
                 (115, -0.7), # right, bottom
                 (0., 0.), # ignored
                 ]
        
        path = Path(verts, codes)
        #add box around area of interest
        ax = plt.subplot()
        patch = patches.PathPatch(path, facecolor='#CCFFCC', lw=1)
        ax.add_patch(patch)
        plt.axhline(y=0.0,xmin=0,xmax=3,c="black",linewidth=1.0,zorder=0)
        plt.savefig(str(pca_modes)+"_ANM_PCA_DIFF.pdf", bbox_extra_artists=(lgd, ), bbox_inches='tight')
        
    def print_from_file(self, modesfile):
        print("Printing modes from file", modesfile)
        #print(len(modes)
        import matplotlib.pyplot as plt
        from matplotlib.path import Path
        import matplotlib.patches as patches
        modes, atm_group = prody.parseNMD(modesfile)

        title=modesfile.split(".")
        modes.setTitle(title[0])
        print_offset = 1
        #print_offset = 738
        plt.rcParams.update({'font.size': 24}) 
        plt.rcParams['xtick.major.pad'] = 10
        plt.rcParams['ytick.major.pad'] = 10
        
        for count in range(6):
            #array with lenghtes of vectors
            a3d = (modes[count].getArrayNx3()**2).sum(axis=1)**0.5
            show = plt.plot(a3d[:], label=("Mode "+str(count+1)))
        plt.xlabel('Residue index')
        plt.ylabel('Lenght of fluctuation vector')
        
        locs,labels = plt.xticks()
        new_labels = ['%d' % (a+print_offset) for a in locs]
        plt.xticks(locs,new_labels)
        plt.ylim([0,0.8])
        plt.xlim([print_offset, modes.numAtoms()])
        
        lgd = plt.legend(loc='center right', bbox_to_anchor=(1.6,0.5))
        plt.savefig(str(modes.getTitle())+"_mode_1_to_6.pdf", bbox_extra_artists=(lgd, ), bbox_inches='tight')
        
        print("Printing modes to file finished")

def parse_Arguments():
    '''
    parsing of command line arguments
    '''
    parser = argparse.ArgumentParser(
        description='Calculate PCA and ANM for given structures',
        epilog="", formatter_class=SmartFormatter)
    parser.add_argument("--pdb", help="""R|specify the pdbs which should be used. 
            Important!!!: They will be overwritten by Prody to match its PDB format,
            in particular histidine names.
            If not specified, all PDBs in the current folder are taken by default.
            Possibilities:
            use all PDBs in directory: --pdb all
            give a list of PDBs, either single or trajectory file: --pdb "name1.pdb name2.pdb name3.pdb name4.pdb"
            """, default='all')
    parser.add_argument("--selection", help="""R|Specify the selection for normal mode analysis
            valid options:
            --selection "calpha": 
            --selection "backbone"
            --selection "all"
    """, default='calpha')
    parser.add_argument("--compare", help="compare PCA and ANM", action="store_true", default=False)
    parser.add_argument("--vmd", help="Visualize result in VMD", action="store_true", default=False)
    parser.add_argument("--debug", help="Print debug information", action="store_true", default=False)
    parser.add_argument("--ref", help="Select a reference structure for PCA and ANM", default='none')
    parser.add_argument("--print_mode", help="Print the square fluctuations of the modes", action="store_true", default=False)
    parser.add_argument("--print_diff", help="Print difference between ANM and PCA modes", action="store_true", default=False)
    parser.add_argument("--print_file", help="Print modes from a NMD file. This will exit the program after printing", default='none')
    parser.add_argument("--anm_only", help="Only perform ANM analysis", action="store_true", default=False)
    return parser.parse_args()

class SmartFormatter(argparse.HelpFormatter):
    '''
    Helper class for argparse
    help beginning with R| will be parsed with new lines
    '''
    def _split_lines(self, text, width):
        # this is the RawTextHelpFormatter._split_lines
        if text.startswith('R|'):
            return text[2:].splitlines()
        return argparse.HelpFormatter._split_lines(self, text, width)


def main(pdb, selection='calpha', vmd=False, compare=False, ref="none", debug=False):

    print("\n\nWarning!!! Reference PDB file, if provided, will be rewritten in"
          " ProDy PDB format, with HID and HIE residues renamed to HIS.\n\n")
    pdb = " ".join(pdb)
        
    #create PCA_analysis object
    Analysis = PCA_ANM_Analysis()
    Analysis.set_default_based_on_argparse(selection, vmd, compare, ref, debug, pdb)
    
    pdbs = Analysis.getPDBs()
    ensemble = Analysis.createEnsemble(pdbs)        
    #set weights for the PCA analysis
    #ensemble = Analysis.setWeights(ensemble)
    ensemble = Analysis.setSelections(ensemble)
    
    pca, pca_modes = Analysis.calcPCA(ensemble)
    anm = Analysis.calcANM(ensemble.getConformation(0))
    
    if Analysis.compare:
        Analysis.compare_modes(anm,pca)
   
    return pca_modes
    
    
if __name__ == '__main__':
    #parse command line options and create help
    args = parse_Arguments()
    main(args.pdb, args.selection, args.vmd, args.compare, args.ref, args.debug)

