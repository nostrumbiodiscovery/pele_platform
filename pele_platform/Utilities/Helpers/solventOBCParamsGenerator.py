'''
Created on Oct 28, 2013

@author: israel

Script to generate OBC parameters for Heteroatoms. The parameters have been extracted from Tinker Molecular package.
This script uses the impact template generated for PELE to identify the atom type and assign the right parameters for the OBC solvent.
If one parameter have not been defined in the OBC templates it puts the default parameter and prints a warning.

Usage:

python solventOBCParamsGenerator.py ALA  

Results:
The ouput will be ALA_OBCParams.txt in the same folder of ALA. Then, to put these new parameters in PELE you have to copy
the data inside the file in the pele template. (Data/OBC/solventParamsHCTOBC.txt) and that's all.

'''

import sys


paramtersLst =[ ['CW', '1.875', '0.72'] ,
                ['NC', '1.7063', '0.79'] ,
                ['CM', '1.875', '0.72'] ,
                ['C*', '1.875', '0.72'] ,
                ['H1', '1.25', '0.85'] ,
                ['CT', '1.9', '0.72'] ,
                ['N2', '1.7063', '0.79'] ,
                ['N*', '1.7063', '0.79'] ,
                ['CR', '1.875', '0.72'] ,
                ['HO', '1.05', '0.85'] ,
                ['NB', '1.7063', '0.79'] ,
                ['H2', '1.25', '0.85'] ,
                ['S', '1.775', '0.96'] ,
                ['NA', '1.7063', '0.79'] ,
                ['H4', '1.25', '0.85'] ,
                ['HC', '1.25', '0.85'] ,
                ['C', '1.875', '0.72'] ,
                ['OH', '1.535', '0.85'] ,
                ['CQ', '1.875', '0.72'] ,
                ['CK', '1.875', '0.72'] ,
                ['O2', '1.48', '0.85'] ,
                ['OS', '1.535', '0.85'] ,
                ['SH', '1.775', '0.96'] ,
                ['HA', '1.25', '0.85'] ,
                ['CB', '1.875', '0.72'] ,
                ['H5', '1.25', '0.85'] ,
                ['CN', '1.875', '0.72'] ,
                ['P', '1.87', '0.86'] ,
                ['N3', '1.625', '0.79'] ,
                ['HP', '1.25', '0.85'] ,
                ['N', '1.7063', '0.79'] ,
                ['H', '1.15', '0.85'] ,
                ['HS', '1.25', '0.85'] ,
                ['CV', '1.875', '0.72'] ,
                ['CA', '1.875', '0.72'] ,
                ['O', '1.48', '0.85'] ,
                ['CC', '1.875', '0.72'] ,]
              #  ['HWS', '1.05', '0.85'] ,
              #  ['OWS', '1.535', '0.85'] ,]

# Parameter List extracted from TINKER source code for OBC  (ksolv.f: +/-line 423)
atomTypesOverlapFactors = [['H','1.25'],
                           ['Li', '1.432'] ,
                           ['C','1.90'],
                           ['N','1.7063'],
                           ['O','1.535'],
                    ['F', '1.47'] ,
                    ['FE', '2.00'] ,   #default parameters
                    ['Ne', '1.39'] ,
                    ['Na', '1.992'] ,
                    ['Mg', '1.70'] ,
                    ['Si', '1.80'] ,
                    ['P', '1.87'] ,
                    ['S', '1.775'] ,
                    ['Cl', '1.735'] ,
                    ['Ar', '1.70'] ,
                    ['K', '2.123'] ,
                    ['Ca', '1.817'] ,
                    ['Br', '1.90'] ,
                    ['Kr', '1.812'] ,
                    ['Rb', '2.26'] ,
                    ['I', '2.10'] ,
                    ['Xe', '1.967'] ,
                    ['Cs', '2.507'] ,
                    ['Ba', '2.188'] ,
                    ['Pt', '2.0'] ,
                  ]

# Parameter List extracted from TINKER source code for OBC  (ksolv.f: +/-line 423)
atomTypesHCTradii = [['H','0.85'],
                     ['C','0.72'],
                     ['N','0.79'],
                     ['O','0.85'],
                     ['F','0.88'],
                  ['P','0.86'],
                  ['S','0.96'],
                  ['Pt','0.80'],
                  ['FE','0.88'],
                  ]


# Parse the impact template to get the residue name and the atom names and atom types
def parseImpactTemplate(impactTemplate):
    
    try:
        fileout = open(impactTemplate)
    except OSError:
        env.logger.info('Impossible to open the file ....  ',impactTemplate)
        sys.exit()
    
    name = ''
    atomNamesAndTypes = []
    atomNumber = []    

    read = False
    bonds = []

    for line in fileout:
        if line[0]!=' ' and len(line.split())==6:
            name = line.split()[0].upper()
            
        if len(line)== 68 and len(line.split())==9:
            tmp = line.split()
            atomNamesAndTypes.append([tmp[4].replace('_',''),tmp[3]])
            atomNumber.append(tmp[0])           

        if 'BOND' in line[:5]: 
            read = True
            continue    
        if 'THET' in line[:5]:
            read = False    
            continue
        if read:
            tmp2 = line.split()
            bonds.append([tmp2[0],tmp2[1]]) 

    # get number of atoms connected
 
    numberOfConnections = []


    for ele in atomNumber:
        counter = 0
        atomAttached = ''
        for bond in bonds:
            if ele==bond[0]:
                counter +=1
                if getShortName(atomNamesAndTypes[int(ele[0])-1][0])[0]=='H':
                    atomAttached= getShortName(atomNamesAndTypes[int(bond[1])-1][0])[0]
            if ele==bond[1]:
                counter +=1
                if getShortName(atomNamesAndTypes[int(ele[0])-1][0])[0]=='H':
                    atomAttached= getShortName(atomNamesAndTypes[int(bond[0])-1][0])[0]

        numberOfConnections.append([counter,atomAttached])

 
    return name, atomNamesAndTypes,numberOfConnections


# Assign overlapFactors and HCT radii using the atom name
def getOverlapscalefactorsFromAtomName(atomName,atomTypesOverlapscalefactors,atomTypesHCTradiinum,numberOfBonds):
   
    shortName = getShortName(atomName)

    radii ='0.80'
    found = False

    overlapFactor, found = assignOverlapFactor(shortName,atomTypesOverlapscalefactors,numberOfBonds,atomName)

    if not found:
        shortName = shortName[0]

        overlapFactor,found = assignOverlapFactor(shortName,atomTypesOverlapscalefactors,numberOfBonds,atomName)

    for atomType in atomTypesHCTradii:
        if shortName==atomType[0].upper():
            radii = atomType[1]

    return overlapFactor,radii,found

def assignOverlapFactor(name,atomTypeOverlapFactorTable,numberOfBonds,realName):
    
    overlapFactor = ''

    found = False

    for atomType in atomTypeOverlapFactorTable:
            if name==atomType[0].upper():
                overlapFactor = atomType[1]
                found = True

    if name=='H' and numberOfBonds[1]=='O': 
        overlapFactor='1.05' 
    if name=='H' and numberOfBonds[1]=='N': overlapFactor='1.15' 

    if name=='C' and numberOfBonds[0]==3: overlapFactor='1.875' 
    if name=='C' and numberOfBonds[0]==2: overlapFactor='1.825' 
    if name=='N' and numberOfBonds[0]==4: overlapFactor='1.625' 
    if name=='N' and numberOfBonds[0]==1: overlapFactor='1.60' 
    if name=='O' and numberOfBonds[0]==1: overlapFactor='1.48' 
 
    return overlapFactor,found

# extract atom type from atom name for hetero atoms
def getShortName(name):
#    shortName = ''.join(i for i in name if not i.isdigit()).upper()

    size = len(name)
    start = 0
    end = 0

    for iterator in range(0,size):
        char = name[iterator]
        if not char.isdigit():
            start = iterator
            break   

    for iterator in reversed(range(0,size)):
        char = name[iterator]
        if not char.isdigit():
            end = iterator
            break   

    short = name[start:end+1]

    return short


#generate the final output template
def generateSolventTemplate(residueName, atomNamesAndTypes,paramtersLst,atomTypesOverlapFactors,atomTypesHCTradii,numberOfBonds):
    
    found = False
    foundNoStd = False
    
    templateParameters = []
    
    for i,ele in enumerate(atomNamesAndTypes):
        found = False
        for param in paramtersLst:
            if ele[1]==param[0]:
                templateParameters.append(residueName+'Z   '+ele[0]+'   '+ele[1]+'    '+param[1]+'   '+param[2]+'\n') 
                found = True
                
        if not found:
            foundNoStd = False
            
            overlapFactor, radii, foundNoStd = getOverlapscalefactorsFromAtomName(ele[0],atomTypesOverlapFactors,atomTypesHCTradii,numberOfBonds[i])
            
            if not foundNoStd:
                env.logger.info('Parameter NOT found in the template database ....... '+ele[0]+'   '+ele[1]+'  using default parameters')

            templateParameters.append(residueName+'Z   '+ele[0]+'   '+ele[1]+'    '+overlapFactor+'   '+radii+'\n') 
        
    return templateParameters

def main(template, out_file):
    # Parameter List extracted from PELE solvent templates for OBC 
    residueName, atomNamesAndTypes, numberOfBonds = parseImpactTemplate(template)
    
    templateParameters =generateSolventTemplate(residueName, atomNamesAndTypes,paramtersLst,atomTypesOverlapFactors,atomTypesHCTradii,numberOfBonds)
    
    fileout = open(out_file, 'a')
    fileout.writelines([item for item in templateParameters]) 
    fileout.close()
    
    

