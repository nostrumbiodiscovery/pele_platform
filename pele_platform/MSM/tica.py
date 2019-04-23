from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import range
from io import open
import matplotlib.pyplot as plt
import os
import glob
import pickle
import numpy as np
import argparse
import itertools
import pyemma.msm as msm
import pyemma.plots as mplt
from pele_platform.Utilities.AdaptivePELE.freeEnergies import cluster, extractCoords
from pele_platform.Utilities.AdaptivePELE.utilities import utilities
from pele_platform.Utilities.AdaptivePELE.atomset import atomset
import pyemma.coordinates as coor
import matplotlib.pyplot as plt
try:
    plt.style.use('ggplot')
except NameError:
    pass

def parse_arguments():
    """
        Create command-line interface
    """
    desc = "Estimate the time-structure based Independent Components (TICA) from a simulation"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("nTICs", type=int, help="Number of Indepent Components to consider")
    parser.add_argument("numClusters", type=int, help="Number of clusters to create")
    parser.add_argument("ligand_resname", type=str, help="Name of the ligand in the PDB")
    parser.add_argument("lag", type=int, help="Lagtime to use in the TICA model")
    parser.add_argument("nTraj", type=int, help="Number of real trajectories per epoch (i.e number of processors-1)")
    parser.add_argument("totalSteps", type=int, default=0, help="Total number of steps in traj. Equivalent to epoch length in adaptive runs")
    parser.add_argument("-o", default=None, help="Path of the folders where the trajectories are stored")
    parser.add_argument("-stride", type=int, default=1, help="Stride, e.g. select one conformation out of every x, default 1, that is take all")
    parser.add_argument("-atomId", type=str, default="", help="Atoms to use for the coordinates of the conformation, if not specified use the center of mass")
    parser.add_argument("-r", "-repeat", action="store_true", help="Force the extraction and repeat the coordinates")
    parser.add_argument("--plot", action="store_true", help="Wheter to make plots")
    parser.add_argument("-top", type=str, default=None, help="Topology file for non-pdb trajectories")
    args = parser.parse_args()
    return args.nTICs, args.numClusters, args.ligand_resname, args.lag, args.nTraj, args.totalSteps, args.o, args.stride, args.atomId, args.r, args.plot, args.top


def get_coords(conformation, atom, lig_name):
    pdb_obj = atomset.PDB()
    if atom:
        atom_name = atom.split(":")[1]
        pdb_obj.initialise(conformation, atomname=atom_name)
        # getAtomCoords returns an array, while getCOM returns a list
        return pdb_obj.getAtom(atom).getAtomCoords().tolist()
    else:
        pdb_obj.initialise(conformation, resname=lig_name)
        return pdb_obj.getCOM()


def find_representative_strucutures(folders, numClusters, nTraj, clusterCenters, projectedUniq, dtrajs):
    centersInfo = {x: {"structure": None, "minDist": 1e6} for x in list(range(numClusters))}
    for i, epoch in enumerate(folders):
        for iTraj, traj in enumerate(projectedUniq[i*nTraj:(i+1)*nTraj]):
            for nSnap, snapshot in enumerate(traj):
                clusterInd = dtrajs[i*nTraj+iTraj][nSnap] 
                dist = np.sqrt(np.sum((clusterCenters[clusterInd]-snapshot)**2))
                if dist < centersInfo[clusterInd]['minDist']:
                    centersInfo[clusterInd]['minDist'] = dist
                    centersInfo[clusterInd]['structure'] = (epoch, iTraj+1, nSnap)
    return centersInfo


def projectTICATrajs(folders, folderPath, ligand_resname, atomId, stride_conformations, nTICs, tica, writeFiles=True, topology=None):
    if writeFiles:
        utilities.makeFolder("tica_COM")
    trajsUniq = []
    projectedUniq = []
    for epoch in folders:
        trajFiles = glob.glob(os.path.join(folderPath, "%s/extractedCoordinates/coord*" % epoch))
        trajFiles.sort(key=lambda x: int(x[x.rfind("_")+1:-4]))
        for trajName in trajFiles:
            trajNum = int(trajName[trajName.rfind("_")+1:-4])
            trajFile = glob.glob(os.path.join(folderPath, "%s/trajectory_%d.pdb" % (epoch, trajNum)))[0]
            snapshotsPDB = utilities.getSnapshots(trajFile, topology=topology)
            trajCOM = [get_coords(snapshot, atomId, ligand_resname) for snapshot in itertools.islice(snapshotsPDB, 0, None, stride_conformations)]
            trajsUniq.append(trajCOM)
            trajLoad = np.loadtxt(trajName)
            if len(trajLoad.shape) == 1:
                trajLoad = trajLoad[np.newaxis, :]
            #De totes agafa les 3 primeres columnes totes les files
            projectedTraj = tica.transform(trajLoad[::stride_conformations])[:, :nTICs]
            projectedUniq.append(projectedTraj)
            if writeFiles:
               np.savetxt("tica_COM/traj_%s_%d.dat" % (epoch, trajNum), np.hstack((np.array(trajCOM), projectedTraj)),
                           header="COM coordinates x\ty\tz\t TICA coordinates\t"+"\t".join(["TICA %d" % tic for tic in range(nTICs)]) + "\n")
    return trajsUniq, projectedUniq


def writeCentersInfo(centersInfo, folderPath, ligand_resname, nTICs, numClusters, trajsUniq, clustersCentersFolder, nTraj, clustersFiles, topology=None):
    if topology is not None:
        topology_contents = utilities.getTopologyFile(topology)
    else:
        topology_contents = None
    if not os.path.exists(clustersCentersFolder):
        os.makedirs(clustersCentersFolder)
    COM_list = []
    for clusterNum in centersInfo:
        epoch, trajNum, snap = centersInfo[clusterNum]['structure']
        COM_list.append(trajsUniq[int(epoch)*nTraj+(trajNum-1)][snap])
        trajFile = glob.glob(os.path.join(folderPath, "%s/trajectory_%d.pdb" % (epoch, trajNum)))
        trajFile = trajFile[0]
        snapshots = utilities.getSnapshots(trajFile, topology=topology)
        pdb_object = atomset.PDB()
        pdb_object.initialise(snapshots[snap], resname=ligand_resname)
        pdb_object.writePDB(str(os.path.join(str(clustersCentersFolder), clustersFiles.format(clusterNum))))

    distances = [[nC, centersInfo[nC]['minDist']] for nC in range(numClusters)]
    np.savetxt(os.path.join(clustersCentersFolder, "clusterDistances_%dcl_%dTICs.dat" % (numClusters, nTICs)), distances)
    utilities.write_PDB_clusters(COM_list, os.path.join(clustersCentersFolder, "clusters_%d_KMeans_allSnapshots.pdb" % (numClusters)))


def make_TICA_plot(nTICs, projected):
    plt.rcParams.update({'legend.markerscale': 10})
    # coords = np.array(projected)
    states = list(range(nTICs))
    for state in states:
        plt.figure()
        # plt.plot(coords[:,:,state].flatten(), 'x', markersize=0.5, label="Tica %d" % (state+1))
        plotNum = 0
        for traj in projected:
            try:
                plt.plot(list(range(plotNum, plotNum+traj.shape[0])), traj[:, state], 'x', markersize=0.5, color="r")
                plotNum += traj.shape[0]
                # plt.plot(traj[:, 2], traj[:, state], 'x', markersize=0.5, color="r")
            except IndexError:
                plt.plot([plotNum], traj[state], 'x', markersize=0.5, color="r")
                plotNum += 1
                # plt.plot(traj[2], traj[state], 'x', markersize=0.5, color="r")
        plt.title("Tica %d" % (state+1))
        # plt.title("Comparing different Tica")
        # plt.xlabel("Tica 3")
        # plt.ylabel("Tica %d" % (state+1))
        # plt.savefig("tica_3_%d_IC.png" % (state + 1))
        plt.savefig("tica_%d_IC.png" % (state + 1))
    # plt.show()
    # import sys
    # sys.exit()


def make_TICA_decomposition(ticaObject, folders, folderPath, lag, overWriteObject=False, kinetic_map=True, commute_map=False):
    if overWriteObject or not os.path.exists(ticaObject):
        trajs = []
        for epoch in folders:
            trajFiles = glob.glob(os.path.join(folderPath, "%s/repeatedExtractedCoordinates/coord*" % epoch))
            trajFiles.sort(key=lambda x: int(x[x.rfind("_")+1:-4]))
            for traj in trajFiles:
                trajs.append(np.loadtxt(traj))
        tica = coor.tica(data=trajs, lag=lag, kinetic_map=kinetic_map, commute_map=commute_map)
        with open(ticaObject, "wb") as f:
            pickle.dump(tica, f)
    else:
        with open(ticaObject, "rb") as f:
            tica = pickle.load(f)
    return tica


def write_TICA_trajs(trajectoryFolder, projected, trajectoryBasename, folders, nTraj):
    if not os.path.exists(trajectoryFolder):
        os.makedirs(trajectoryFolder)

    for i, epoch in enumerate(folders):
        for iTraj, traj in enumerate(projected[i*nTraj:(i+1)*nTraj]):
            auxArr = np.zeros_like(traj[:, 0])
            # Add a first column of indexes because it is the format that the
            # cluster module of the freeEnergie module reads
            np.savetxt(os.path.join(trajectoryFolder, "%s_%d_%d.dat" % (trajectoryBasename[:-1], int(epoch), iTraj+1)), np.hstack((auxArr.reshape(-1, 1), traj)))


def cluster_TICA_space(numClusters, trajectoryFolder, trajectoryBasename, stride, clusterCountsThreshold):
    clusteringObject = cluster.Cluster(numClusters, trajectoryFolder,
                                       trajectoryBasename, alwaysCluster=False,
                                       stride=stride)
    clusteringObject.clusterTrajectories()
    clusteringObject.eliminateLowPopulatedClusters(clusterCountsThreshold)
    return clusteringObject


def main(nTICs, numClusters, ligand_resname, lag, nTraj, n_steps, out_path=None, stride_conformations=1, atomId="", repeat=False, plotTICA=False, topology=None,
    clustersCentersFolder="clustersCenters", trajectoryFolder="tica_projected_trajs", trajectoryBasename="tica_traj*", stride=1, clusterCountsThreshold=0, 
    ticaObject = "tica.pkl", clustersFiles="clusters_{}.pdb"):

    lags = [1,2,5,10,20,50,100,200]
    if out_path is None:
        folderPath = ""
        curr_folder = "."
    else:
        folderPath = out_path
        curr_folder = out_path

    folders = utilities.get_epoch_folders(curr_folder)
    folders.sort(key=int)

    if not os.path.exists(os.path.join(folderPath, "0/repeatedExtractedCoordinates/"))or repeat:
        # Extract ligand and alpha carbons coordinates
        extractCoords.main(folder_name=curr_folder, lig_resname=ligand_resname, numtotalSteps=n_steps, protein_CA=False, non_Repeat=False, sidechains=True, sidechain_folder="../output_clustering/initial*", nProcessors=20, parallelize=True)

    tica = make_TICA_decomposition(ticaObject, folders, folderPath, lag)

    # Select the desired number of independent components from the full
    # decomposition
    Y = tica.get_output(list(range(nTICs)))
    print('number of trajectories = ', np.shape(Y)[0])
    print('number of frames = ', np.shape(Y)[1])
    print('number of dimensions = ',np.shape(Y)[2])
    write_TICA_trajs(trajectoryFolder, Y, trajectoryBasename, folders, nTraj)
    clustering = coor.cluster_kmeans(Y, k=250, max_iter=1000)
    dtrajs = clustering.dtrajs
    mplt.plot_free_energy(np.vstack(Y)[:,0], np.vstack(Y)[:,1])
    cc_x = clustering.clustercenters[:,0]
    cc_y = clustering.clustercenters[:,1]
    plt.plot(cc_x,cc_y, linewidth=0, marker='o', markersize=5, color='black')
    plt.show()
    its = msm.timescales_msm(dtrajs, lags=800, nits=800)
    mplt.plot_implied_timescales(its, ylog=False, units='steps', linewidth=2)
    plt.show()
    #clusteringObject = cluster_TICA_space(numClusters, trajectoryFolder, trajectoryBasename, stride, clusterCountsThreshold)
    #trajsUniq, projectedUniq = projectTICATrajs(folders, folderPath, ligand_resname, atomId, stride_conformations, nTICs, tica, topology=topology)
    #its = msm.its(projectedUniq, lags=lags)
    #clusterCenters = clusteringObject.clusterCenters
    #dtrajs = clusteringObject.assignNewTrajectories(projectedUniq)
    #centersInfo = find_representative_strucutures(folders, numClusters, nTraj, clusterCenters, projectedUniq, dtrajs)
    #writeCentersInfo(centersInfo, folderPath, ligand_resname, nTICs, numClusters, trajsUniq, clustersCentersFolder, nTraj, clustersFiles, topology=topology)
    #if plotTICA:
    #    make_TICA_plot(nTICs, projected)
    #return clusteringObject

if __name__ == "__main__":
    nICs, n_clusters, lig_resname, lagtime, n_traj, n_step, output_path, stride_conformation, atomIds, repeat_true, plots, top = parse_arguments()
    main(nICs, n_clusters, lig_resname, lagtime, n_traj, n_step, output_path, stride_conformation, atomIds, repeat_true, plots, topology=top)
