"""
$Revision: 1.13 $
For ligands:
Reads in a maestro mae file and makes a "rotamer enabled" template and the rotamer library to accompany it.  
This consists of finding the backbone core that results in the least number of child bonds rotated with any rotatable bond rotation. 
Reads in the rotatable bonds from a macromodel atomtyping (should be easily modifyiable to read them from stdin).  
Hetgrp_ffgen is used for atomtyping and to determine the geometry in the template fromt the mae.  The mae must only have the residue to atomtype in it. 

Builds rotamer libraries for arbitrary ligand molecules by building and combining rotamer libraries.  There are two general algorithms that are implemented.  This first is using macromodel (or some other external tool) to sample the degrees of freedom and converting the resulting ensemble into a rotamer library.  The second is identifying the rotatable bonds, again using macromodel, and assigning rotamer libraries to these bonds.  For most bonds a simple freely rotatable library ( 0,10,20,30...360), but in the case of closed rings special libraries are built using macromodel sampling.  These component rotamer libraries are then arranged into groups for use in PLOP.  Each group consists of a single tree rooted at the central core.  This core can either be used chosen, or will be chosen based on an algorithm that minimizes the number of bond lengths from the farthest leeaf to the trunk.  Any built rotamer libraries are placed in the current directory and a file called <resname>.rot.assign is also written.  This tells PLOP how to assemble the full combinatoric library that will be used in sidehchain prediction/monte carlo.     

For unnatural amino acids:
Requires: 1) a maestro mae file of an unnatural amino acid with no NMA or ACE termini; the N-H and C=0 are to be left as they are found in a peptide
          2) a template file (specified by -t=<FILENAME>) created by hetgrp-ffgen using a maestro mae file of an unnatural amino acid with NMA or ACE termini present
Option Settings Required:  -unnat=yes -t=<FILENAME> [other options] <RESIDUE_WITHOUT_CAPPING_GROUPS>.mae
Outputs: 1) a re-ordered template file (old one is backed up in FILENMAE.hetgrp_ffgen)
         2) a PLOP nonstandard residue specification for pasting into a PLOP control file, both to stdout and to <maefile>_torsions.txt


Options:

   --core <an1>      Give one atom of the core section

   --n <number>      Maximum Number of Entries in Rotamer File

   --mtor <number>   Gives the maximum number of torsions allowed in each
                     group.  Will freeze bonds to extend the core if 
                     necessary

   --clean           Clean Intermiadiate files



Mae file should be properly atomtyped

Most common problem:  As part of this procedure the pdb atom names are often renamed when they are not unique. Alsothis procedure works best if the ligand is minimized.  For these two reasons an atomtyped, minimzed version of the input ligand is written to (input).PlopRotTemp.pdb.  If at all possible, use the ligand structure and atom names from this file in any subsequent plop runs.    

examples:
Build a rotamer library for the following ligand at a grid resolution of 20 deg using PLOP/PRIME to combine libraries
$SCHRODINGER/utilities/python PlopRotTemp.py 3ert_lig.mae -g=20

Build a rotamer library for the following ligand at using CGEN sampling in macromodel.
$SCHRODINGER/utilities/python PlopRotTemp.py 3ert_lig.mae -a=CGEN

Build a rotamer library for the following ligand at a grid resolution of 20 deg using macromodel to sample any rings and combining this with freely rotatable libraries in PLOP to create combined libraries for the ligand.
$SCHRODINGER/utilities/python PlopRotTemp.py 1rth_lig.mae -r=yes

Make libraries for rotatable bonds in ligand.mae up to a maximum of 4 rotatable bonds in each library
All additional bonds are defined as backbone and are sampled with CGEN to produce a backbone libary
$SCHRODINGER/utilities/python PlopRotTemp.py ligand.mae -mtor=4 -ba=CGEN

For a given ligand named LIG the following files will be created:
lig                - Template file for use in PLOP, its zmatrix matches the libraries created
LIG.rot.assign     - Summary of all libraries build or used for this ligand read into plop with the command
                     "rot assign all"
LIG???.side        - (OPTIONAL) Component sidechains libraries created if there are closed rings or CGEN sampling is used
LIG__B.back        - (OPTIONAL) Backbone sidechain library

---------------------------------------------------------------------------------------------------------



All jobs run on the localhost


USAGE: "$SCHRODINGER/utilities/python main.py [file.mae] --options valueoption"

HELP: $SCHRODINGER/utilities/python main.py --help

"""


import argparse
import sys
import os
import re
import PlopRotTemp as pl
import schrodinger.application.macromodel.utils as mu
import schrodinger.application.macromodel.tools as mt
import schrodinger.job.jobcontrol as jc
import schrodinger.infra.mm as mm




#Defaults
template_file = ""
debug = 0  # 1 means don't run exteral commands (assumes output is already there)
conf_file = 'none';
output_template_file = ""
gridres = "10.0"
nsamp = 10000
nrot = 1000
Ecut = 100
use_rings = 0
do_init_min = 1
user_core_atom = -1
max_dist_eq = 0.25
user_tors = []
back_tors = []
back_algorithm = "none"
back_conf_file = ""
hetgrp_opt = ""
use_mae_charges = 0
OPLS = "2005"
max_tors = 5
user_fixed_bonds = []
files2clean = []
use_mult_lib = 1
run_conf = 0
algorithm = "MCMM"
clean = False
#gridres_oh = ""
gridres_oh = gridres
unnat_res = 0  # for old-style PLOP nonstandard side chain
resno = 1  #for old-style PLOP nonstandard side chain
chain = 'A'  #for old-style PLOP nonstandard side chain
grow = 0
tree = 0  # for old-style PLOP nonstandard ligand tree-style torsion reordering
R_group_root_atom_name = 'None'  # which atom do you want to start sampling at?
# R_group_root_atom_name added for "sar" option for only sampling after this atom has been passed
# a hack for post-reordering assignment of rank 0 to reordered atoms between 
#the core and the atom with PDB atom name = R_group_root_atom_name
# only enabled for library build-up in PLOP, not for macromodel sampling
# careful with multiple groups -- will only work properly if these are after the R group atom in the reordered template file



parser = argparse.ArgumentParser()
parser.add_argument("mae_file", type=str, help="ligand maestro mae file")
parser.add_argument("--core", type=int, help="Give one atom of the core section")
parser.add_argument("--mtor", type=int, help="Gives the maximum number of torsions allowed in each \
                              group.  Will freeze bonds to extend the core if \
                              necessary.")
parser.add_argument("--n", type=int, help="Maximum Number of Entries in Rotamer File")
parser.add_argument("--clean", help="Whether to clean up all the intermediate files", action='store_true')
args = parser.parse_args()

if args.mae_file:
  if(args.mae_file.endswith('.mae')):
    mae_file = args.mae_file
  else:
    raise Exception('A .mae file is needed')
if args.core:
  print('\nUsing user core information : {}\n'.format(args.core))
if args.mtor:
  max_tors = args.mtor
  if max_tors>5:
    raise Exception('Maximum mTor number 5')
  print('\nUsing {} as a maximum number of Rotamers\n'.format(max_tors))
if args.n:
  nrot = args.n
  print('\nUsing {} as a Maximum Number of Entries in Rotamer Files\n'.format(nrot))
if args.clean:
  clean = args.clean



#########################COMENT#################################

# Process options
"""
if (gridres_oh == ""): 
    gridres_oh = gridres
if (use_mae_charges == 1):
    hetgrp_opt = hetgrp_opt + '-use_mae_charges'

if (run_conf == 0): 
    conf_file = 'none'
"""
#######################COMENT################################33

####################CHANGE MACROMODEL###########################
if ((user_tors != [] or user_fixed_bonds != []) and conf_file == ''):
    raise Exception("Cannot call Macromodel to perform sampling with user defined torsions")
####################CHANGE MACROMODEL###########################

if (unnat_res == 1):
    init_min = 0  #the input mae file is for a peptide and will not have a suitable Lewis structure
    if (template_file == ""):
        print("Cannot use unnatural residues without pre-made template files!")
        sys.exit(-1)
    use_mult_lib = 1  # so a dummy conformational search is performed just to see which bonds are rotatble
    use_rings = 0  #for now; I'm not sure that ring torsions will follow the tree pattern appropriately, this would need testing
    #For now, just try low energy ring conformations

###############POT COMENTAR#################
if (R_group_root_atom_name != 'None'):
    use_mult_lib = 1  # so a dummy conformational search is performed just to see which bonds are rotatble
    use_rings = 0  #for now; I'm not sure that ring torsions will follow the tree pattern appropriately, this would need testing
    #For now, just try low energy ring conformations
###############POT COMENTAR#################

####################REMOVE MACROMODEL###########################
# Create ComUtil instance , define potential energy surface: solution phase, OPLSAA
# Serial mode enabled so each structure is used to seed a unique search
mcu_conf = mu.ComUtil(ffld='opls2005', serial=True, solv=True, nant=False, demx=True)
mcu_dummy = mu.ComUtil(ffld='opls2005', serial=True, solv=True, nant=False, demx=True)
mxu = mu.CluUtil()


# There are two debug switch associated with AUTO
# Debug output appears in jobname.log
mcu_conf.SOLV[2] = 1  # water
mcu_dummy.DEBG[1] = 520  # Debugging output is read to determine zmatrix
mcu_dummy.DEBG[2] = 521

####################REMOVE MACROMODEL###########################

root = pl.get_root_path(mae_file)

print("\n")
print("INPUT")
print("mae_file {}".format(mae_file))
print("root {}".format(root))
print("OPLS {}".format(OPLS))
print("hetgrp options '{}'".format(hetgrp_opt))
print("User template file '{}'".format(template_file))
print("User output template file '{}'".format(output_template_file))



#########################CHANGE HETGRP_FFGEN ligand preparation for PELE###################
#Build a template file
print("\n")
print("TEMPLATE GENERATION")
[template_file, output_template_file, mae_file_hetgrp_ffgen, files, resname] = \
    pl.build_template(mae_file, root, OPLS, hetgrp_opt, template_file, \
                   output_template_file)
for f in files:
    files2clean.append(f)

#########################CHANGE HETGRP_FFGEN ligand preparation for PELE###################

####################CHANGE MACROMODEL MINIMIZATION OF LIGAND-->PymoChimera, OPenMM###########################
#(.mae)-->pdb (pymol)
#OpenMM (minimize)
#(pymol)-->(.mae)
print("\n")
if (do_init_min == 1):
    mcu_mini = mu.ComUtil(ffld='opls2005', serial=True, solv=True, nant=False, demx=True)
    mcu_mini.SOLV[2] = 1  # water
    mini_root = root + "_mini"
    com_file = mcu_mini.mini(mae_file_hetgrp_ffgen, mini_root + '.com')
    print('\nMINIMIZATION\nRunning minimization: {0} -> {1} -out.mae\n'.format(mae_file_hetgrp_ffgen, mini_root))
    if (not debug):
        cmd = mcu_mini.getLaunchCommand(com_file)
        job = jc.launch_job(cmd)
        job.wait()
        files2clean.append(mini_root + '-out.mae')
        files2clean.append(mini_root + '.log')
#       files2clean.append(mini_root + '-out.tmp')
#        files2clean.append(mini_root + '.com')
    mae_min_file = mini_root + "-out.mae"
else:
    print('\nSkipping Minimization\n ')
    mae_min_file = mae_file_hetgrp_ffgen

#Run the Dummy Conformation Search to Find Bonds
mcu_dummy.MCMM[1] = 1  # Take 1,000 steps
mcu_dummy.MCMM[2] = 1  # Store up to 1000 values
mcu_dummy.MINI[3] = 0  # Don't minimize
mcu_dummy.DEMX[5] = 100
mcu_dummy.DEMX[6] = 500
print('\nCONFORMATIONAL SEARCH\nRunning dummy conformation search to find bonds\n')
com_file = mcu_dummy.mcmm(mae_min_file, root + '_IDbonds.com')
log_file = root + '_IDbonds.log'
if (not debug):
    cmd = mcu_dummy.getLaunchCommand(com_file)
    job = jc.launch_job(cmd)
    job.wait()
    files2clean.append(com_file)
    files2clean.append(log_file)
    files2clean.append(root + '_IDbonds-out.mae')
    files2clean.append(root + '_IDbonds-out.ouL')

#################CHANGE MACROMODEL MINIMIZATION OF LIGAND + CONFORMATIONAL SEARCH###########################


####################SCHRODINGER-->get atoms from rings & see whether or not they are bonded###########################
print("\n")
if (unnat_res == 1):
    [mae_num, parent, rank, tors, use_rings, group, tors_ring_num] = \
        pl.FindCoreAA(mae_min_file, user_fixed_bonds, log_file, use_rings, use_mult_lib, user_core_atom, user_tors)
    tors_ring_num = []
    for t in tors: tors_ring_num.append(0);
else:
    print('FINDING CORE')
    if (grow == 1 and user_core_atom == -1): user_core_atom = -2
    #######Assign_rank--> Extremely slow!
    [mae_num, parent, rank, tors, use_rings, group, back_tors, tors_ring_num] = \
        pl.FindCore(mae_min_file, user_fixed_bonds, log_file, use_rings, \
                 use_mult_lib, user_core_atom, user_tors, back_tors, max_tors, R_group_root_atom_name)
if (use_rings == 1):
    print("Found flexible rings")

newtors = []
if (unnat_res == 1 or grow == 1 ):
    newtors = pl.ReorderTorsionsAA(tors, mae_num)

####################SCHRODINGER###########################



#Change from mae files atom numbering to the template file ones
#Convert Torsions to match the template file atom numbering
#Ring numbers don't have to be changed
[mae2temp, temp2mae] = pl.MatchTempMaeAtoms(mae_min_file, template_file)
old_atom_num = [];
new_tors = [];
new_back_tors = [];
for i in mae_num:
    old_atom_num.append(-100)
for i in range(len(mae2temp)):
    old_atom_num[i] = mae2temp[mae_num[i]]
for i in range(len(tors)):
    temp = [mae2temp[tors[i][0]], mae2temp[tors[i][1]]]
    new_tors.append(temp)
for i in range(len(back_tors)):
    temp = [mae2temp[back_tors[i][0]], mae2temp[back_tors[i][1]]]
    new_back_tors.append(temp)
tors = []
for i in range(len(new_tors)):
    temp = [old_atom_num.index(new_tors[i][0]), old_atom_num.index(new_tors[i][1])]
    tors.append(temp)
back_tors = []
for i in range(len(new_back_tors)):
    temp = [old_atom_num.index(new_back_tors[i][0]), old_atom_num.index(new_back_tors[i][1])]
    back_tors.append(temp)

#Make (or read) original tempalte file
print('\n')
print('CREATE ROTAMER TEMPLATE FILE: {}'.format(output_template_file))
names = pl.ReorderTemplate(old_atom_num, parent, rank, template_file, output_template_file,
                        R_group_root_atom_name=R_group_root_atom_name)

[tors, tors_ring_num, zmat_atoms] = pl.FindTorsAtom(tors, tors_ring_num, parent)
#Eliminate Torsions in the backbone (included when entire rings are included in the torsions)
[tors, tors_ring_num, zmat_atoms] = pl.EliminateBackboneTors(tors, tors_ring_num, zmat_atoms, rank)

if (unnat_res == 1 or grow == 1):
    mynonstandard = pl.TetherRotBonds(mae_file, chain, resno, log_file, newtors)
    mynonstandard.output_rotbonds(R_group_root_atom_name=R_group_root_atom_name)


else:

    # Reorder the torsions
    for i in range(len(tors)):
        tors[i].sort()
    for i in range(len(tors)):
        for j in range(i + 1, len(tors)):
            if (tors[i] > tors[j]):
                temp = tors[i];
                tors[i] = tors[j];
                tors[j] = temp;
                temp = tors_ring_num[i];
                tors_ring_num[i] = tors_ring_num[j];
                tors_ring_num[j] = temp

    [tors, tors_ring_num, zmat_atoms] = pl.FindTorsAtom(tors, tors_ring_num, parent)


################################CHANGE MACROMODEL CONFORMATIONAL SEARCH######################
#Run the conformational Search
conf_root = root + "_conf"
if (conf_file == conf_root + '-out.mae'):
    raise Exception('Must use different name for conformational file')

if (conf_file == ''):
    run_conf = 1
else:
    run_conf = 0
#Can I comment that?
"""
if (run_conf == 1 ):  #We are actually going to run a csearch
    print('Taking {0} steps and storin {1} conformations'.format(nsamp, nrot))
    print('Algorithm to be used is {}'.format(algorithm))
    print('Energy cutoff is {} kcals/mole'.format(Ecut))
    conf_file = conf_root + '-out.mae'
    if (algorithm == "MCMM" or algorithm == "mcmm"):
        mcu_conf.MCMM[1] = nsamp  # Take X steps
        mcu_conf.MCMM[2] = nrot  # Store up to Y values
        mcu_conf.MINI[1] = 1  # PRCG
        mcu_conf.MINI[3] = 50  # iterations of minimze
        mcu_conf.DEMX[5] = Ecut  #cutoffs in kJ/mole
        mcu_conf.DEMX[6] = Ecut * 5  #cutoffs in kJ/mole
        com_file = mcu_conf.mcmm(mae_min_file, conf_root + '.com')
    elif (algorithm == "CGEN" or algorithm == "cgen"):
        mcu_conf.CGEN[1] = nsamp
        mcu_conf.CGEN[2] = nrot
        mcu_conf.MINI[1] = 1  # PRCG
        mcu_conf.MINI[3] = 50  # iterations of minimze
        mcu_conf.DEMX[5] = Ecut  #cutoffs in kJ/mole
        mcu_conf.DEMX[6] = Ecut * 5  #cutoffs in kJ/mole
        com_file = mcu_conf.cgen(mae_min_file, conf_root + '.com')
    else:
        raise Exception("Algorithm ", algorithm, " not recognized\n");
    if (user_tors == []):
        if (run_conf == 1):
            print('Running Conformational Search: {0} -> {1}'.format(mae_min_file, conf_file))
        if (not debug):
            cmd = mcu_conf.getLaunchCommand(com_file)
            job = jc.launch_job(cmd)
            job.wait()
            files2clean.append(conf_root + '.com')
            files2clean.append(conf_root + '.log')
            files2clean.append(conf_root + '-out.mae')
            files2clean.append(conf_root + '-out.ouL')
    else:
        raise Exception("Cannot combine user defined torsions and MM search");

################################CHANGE CONFORMATIONAL SEARCH######################

#Cluster the output
#if(run_conf==1):
#  clust_name=root+'_clust.mae'
#  print "Clustering the results of the conformational search",conf_file,' -> ',clust_name
#  count=mt.count(conf_file)
#  print 'Clustering from ',count,' to ',nrot
#  thresh = count - nrot + 1
#  print 'Found members in cluster',count,nrot
#  if count > nrot and nstore > nrot:
#     mxu.arms[0] = 'heavy'
#     mxu.writerep = " ".join([`thresh`,clust_name,"all"])
#     mxu.doCluster(conf_file)
#  else:
#     os.rename(conf_file,clust_name)
#else:  
#   clust_name=conf_file


################################CHANGE MACROMODEL CONFORMATIONAL SEARCH######################
#Run the conformational Search for the backbone



if (back_tors != [] and back_algorithm != "none"):
    conf_root = root + "_backconf"
    if (back_conf_file == conf_root + '-out.mae'):
        raise Exception('Must use different name for backbone conformational file')

    if (back_conf_file == ''):
        run_conf = 1
    else:
        run_conf = 0
    if (run_conf == 1):  #We are actually going to run a csearch
        back_conf_file = conf_root + '-out.mae'
        print('Running Backbone Conformational Search: ' + mae_min_file + ' -> ' + back_conf_file)
        print('Taking '+ nsamp+ ' steps and storing '+ nrot+ ' conformtations')
        print('Algorithm to be used is {}'.format(back_algorithm))
        print('Energy cutoff is {} kcals/mole'.format(Ecut))
        if (back_algorithm == "MCMM" or back_algorithm == "mcmm"):
            mcu_conf.MCMM[1] = nsamp  # Take X steps
            mcu_conf.MCMM[2] = nrot  # Store up to Y values
            mcu_conf.MINI[1] = 1  # PRCG
            mcu_conf.MINI[3] = 50  # iterations of minimze
            mcu_conf.DEMX[5] = Ecut  #cutoffs in kJ/mole
            mcu_conf.DEMX[6] = Ecut * 5  #cutoffs in kJ/mole
            com_file = mcu_conf.mcmm(mae_min_file, conf_root + '.com')
        elif (back_algorithm == "CGEN" or back_algorithm == "cgen"):
            mcu_conf.CGEN[1] = nsamp
            mcu_conf.CGEN[2] = nrot
            mcu_conf.MINI[1] = 1  # PRCG
            mcu_conf.MINI[3] = 50  # iterations of minimze
            mcu_conf.DEMX[5] = Ecut  #cutoffs in kJ/mole
            mcu_conf.DEMX[6] = Ecut * 5  #cutoffs in kJ/mole
            com_file = mcu_conf.cgen(mae_min_file, conf_root + '.com')
        else:
            raise Exception("Algorithm ", back_algorithm, " not recognized\n");
        if (not debug):
            cmd = mcu_conf.getLaunchCommand(com_file)
            job = jc.launch_job(cmd)
            job.wait()
            files2clean.append(conf_root + '.com')
            files2clean.append(conf_root + '.log')
            files2clean.append(conf_root + '-out.mae')
            if (back_algorithm == "CGEN"):
                files2clean.append(conf_root + '-out.mmo')


    if (unnat_res != 1):
        if (back_conf_file != '' and back_conf_file != 'none'):
            back_lib = resname.upper() + "__B"
            pl.make_lib_from_mae(back_lib, "back", back_conf_file, back_tors, names, \
                              parent, old_atom_num, mae2temp, temp2mae, gridres, gridres)

#if(grow == 1):
#  #Convert conf file to pdbs if necessary
#  pdb_root = root+".PlopRotTemp.pdb"
#  line="$SCHRODINGER/utilities/pdbconvert -imae "+mae_file+" -opdb "+pdb_root+" -num_models 1"
#  file2clean=[]
#  print "Converting mae file to pdb format -> ",pdb_root
#  os.system(line)
################################CHANGE MACROMODEL CONFORMATIONAL SEARCH######################



"""
back_lib = "";
if (unnat_res != 1):  
    if (conf_file != 'none'):
        rotamers_file = pl.make_libraries(resname, conf_file, root, names, zmat_atoms, group, use_rings, use_mult_lib,
                       output_template_file, gridres, debug)
        print("\n")
        print("CREATE ROTAMER LIBRARY")
        print(rotamers_file)
        print("\n")


    ############################CHANGE MACROMODEL--> ring database#########################
    else:
        if (len(zmat_atoms) > 0):
            ring_libs = pl.build_ring_libs(mae_min_file, root, resname, tors, \
                                        tors_ring_num, names, rank, parent, old_atom_num, mae2temp, gridres,
                                        files2clean, debug)
    ############################CHANGE MACROMODEL#########################
        else:
            ring_libs = []
            print("No rotatable sidechains found")
        
        rotamers_file = pl.find_build_lib(resname, mae_min_file, root, tors, names, group, gridres, gridres_oh, use_rings, back_lib,
                              tors_ring_num, ring_libs, debug)
        print("\n")
        print("CREATE ROTAMER LIBRARY")
        print(rotamers_file)
        print("\n")

if (clean):
    print(files2clean)
    for file in files2clean:
        print('Removing Intermediate File: {}'.format(file))
        os.remove(file)
