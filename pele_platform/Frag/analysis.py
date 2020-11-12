from Bio.PDB import PDBParser
import glob, numpy, os
import pandas as pd

parser = PDBParser()


def getAtomFromRef(file, atomname='C29', resname="FMM"):
    reference = parser.get_structure("reference", file)
    for residue in reference.get_residues():
        if residue.resname == resname:
            for atom in residue.get_atoms():
                if atom.name == atomname:
                    atom_ref = atom
                    break
            break
    return atom_ref


def getBestDistance(pdb_file, pointRef, resname="GRW"):
    struct = parser.get_structure("epoch", pdb_file)
    bestAtom = ''
    bestDist = ''
    refVector = numpy.array(pointRef)
    for residue in struct.get_residues():
        if residue.resname == resname:
            for atom in residue.get_atoms():
                atomVector = numpy.array(atom.get_coord())
                dist = numpy.linalg.norm(atomVector - refVector)
                if (bestDist == '' or dist < bestDist):
                    bestAtom = atom
                    bestDist = dist
    return (bestDist, bestAtom)


def get_distance(pointRef, path, resname):
    bestDist = ''
    bestFileDistance = ''
    bestAtom = ''
    for file in glob.glob(path):
        distance, atom = getBestDistance(file, pointRef, resname=resname)
        if (bestDist == '' or distance < bestDist):
            bestDist = distance
            bestFileDistance = file
            bestAtom = atom
    return (bestDist, bestFileDistance, bestAtom)


def get_distance_be(pointRef, path, resname):
    bestDist = ''
    bestFileDistance = ''
    bestAtom = ''
    bestBE = 0
    bestFileBE = ''
    distList = []
    beList = []
    for file in glob.glob(path):
        distance, atom = getBestDistance(file, pointRef, resname=resname)
        distList.append(distance)
        if (bestDist == '' or distance < bestDist):
            bestDist = distance
            bestFileDistance = file
            bestAtom = atom
        be = file.split("_BindingEnergy")[1]
        be = float(be[:-4])
        beList.append(be)
        if (be == 0 or be < bestBE):
            bestBE = be
            bestFileBE = file
    distListNorm, beListNorm = normalize_lists(distList, beList)
    normList = numpy.sqrt(numpy.multiply(numpy.power(distListNorm, 2), numpy.power(beListNorm, 2)))
    bestNormalize = min(normList)
    indexBestNorm = numpy.where(normList == bestNormalize)
    bestNormFile = glob.glob(path)[indexBestNorm[0][0]]
    distBestNorm = distList[indexBestNorm[0][0]]
    beBestNorm = beList[indexBestNorm[0][0]]
    return (
    bestDist, bestFileDistance, bestAtom, bestBE, bestFileBE, bestNormalize, bestNormFile, distBestNorm, beBestNorm)


def normalize_lists(distList, beList):
    maxDist = max(distList)
    distList = [x / maxDist for x in distList]
    minBE = min(beList)
    beList = numpy.power([x / minBE for x in beList], -1)
    return (distList, beList)


def main(refFile='', path="", atomCoords='', resname="GRW", pattern=''):
    if atomCoords == '':
        atomCoords = getAtomFromRef(refFile).get_coord()
    pathSearch = os.path.join(path, "*{}*".format(os.path.splitext(pattern)[0]))
    configFiles = glob.glob(pathSearch)
    files = []
    bestFilesFoundDistance = []
    dists = []
    bestFilesFoundBE = []
    BEs = []
    norm = []
    bestFilesFoundNorm = []
    distNormList = []
    beNormList = []
    for file in configFiles:
        if os.path.isdir(os.path.join(file, "top_result")):
            path = os.path.join(file, "top_result/epoch*")
            bestDist, bestFileDistance, bestAtom, bestBE, bestFileBE, bestNormalize, bestNormFile, distBestNorm, beBestNorm = get_distance_be(
                atomCoords, path, resname)
            distMath = format(bestDist, '.15g')
            files.append(file)
            bestFilesFoundDistance.append(bestFileDistance)
            dists.append(distMath)
            bestFilesFoundBE.append(bestFileBE)
            BEs.append(bestBE)
            norm.append(bestNormalize)
            bestFilesFoundNorm.append(bestNormFile)
            distNormList.append(distBestNorm)
            beNormList.append(beBestNorm)
    dataframe = {'File': files, 'BestFileDistance': bestFilesFoundDistance, 'Distance': dists,
                 'BestFileBE': bestFilesFoundBE, 'BE': BEs, 'BestFileNormalization': bestFilesFoundNorm,
                 'BestNormalization': norm, 'DistanceBestNormalization': distNormList,
                 'BEBestNormalization': beNormList}
    df = pd.DataFrame(dataframe,
                      columns=['File', 'BestFileDistance', 'Distance', 'BestFileBE', 'BE', 'BestFileNormalization',
                               'BestNormalization', 'DistanceBestNormalization', 'BEBestNormalization'])
    df.to_csv('point_analysis.csv', index=False, header=True)


