import sys
import os
import glob
import numpy as np
import argparse
from pele_platform.AdaptivePELE.utilities import utilities
from pele_platform.AdaptivePELE.freeEnergies import cluster, extractCoords
from pele_platform.AdaptivePELE.analysis import splitTrajectory

def parseArgs():
    parser = argparse.ArgumentParser(description="Script that reclusters the Adaptive clusters")
    parser.add_argument('nClusters', type=int)
    parser.add_argument("ligand_resname", type=str, help="Name of the ligand in the PDB")
    parser.add_argument("-atomId", nargs="*", default="", help="Atoms to use for the coordinates of the conformation, if not specified use the center of mass")
    parser.add_argument('-o', type=str, help="Output folder", default="None")
    args = parser.parse_args()
    return args.nClusters, args.ligand_resname, args.atomId, args.o


def writePDB(pmf_xyzg, title="clusters.pdb"):
    templateLine = "HETATM%s  H%sCLT L 502    %s%s%s  0.75%s           H\n"

    content = ""
    for j, line in enumerate(pmf_xyzg):
        number = str(j).rjust(5)
        number3 = str(j).ljust(3)
        x = ("%.3f" % line[0]).rjust(8)
        y = ("%.3f" % line[1]).rjust(8)
        z = ("%.3f" % line[2]).rjust(8)
        g = 0
        content += templateLine % (number, number3, x, y, z, g)

    with open(title, 'w') as f:
        f.write(content)


def writeInitialStructures(centers_info, filename_template, topology=None):
    for cluster_num in centers_info:
        epoch_num, traj_num, snap_num = map(int, centers_info[cluster_num]['structure'])
        trajectory = "%d/trajectory_%d.xtc" % (epoch_num, traj_num) if topology else "%d/trajectory_%d.pdb" % (epoch_num, traj_num)
        snapshots = utilities.getSnapshots(trajectory, topology=topology)
        if not topology:
            with open(filename_template % cluster_num, "w") as fw:
                fw.write(snapshots[snap_num])
        else:
            splitTrajectory.main("", [trajectory, ], topology, [snap_num+1,],template=filename_template % cluster_num)


def get_centers_info(trajectoryFolder, trajectoryBasename, num_clusters, clusterCenters):
    centersInfo = {x: {"structure": None, "minDist": 1e6, "center": None} for x in xrange(num_clusters)}

    trajFiles = glob.glob(os.path.join(trajectoryFolder, trajectoryBasename))
    for traj in trajFiles:
        _, epoch, iTraj = os.path.splitext(traj)[0].split("_", 3)
        trajCoords = np.loadtxt(traj)
        if len(trajCoords.shape) < 2:
            trajCoords = [trajCoords]
        for snapshot in trajCoords:
            nSnap = snapshot[0]
            snapshotCoords = snapshot[1:]
            dist = np.sqrt(np.sum((clusterCenters-snapshotCoords)**2, axis=1))
            for clusterInd in xrange(num_clusters):
                if dist[clusterInd] < centersInfo[clusterInd]['minDist']:
                    centersInfo[clusterInd]['minDist'] = dist[clusterInd]
                    centersInfo[clusterInd]['structure'] = (epoch, int(iTraj), nSnap)
                    centersInfo[clusterInd]['center'] = snapshotCoords
    return centersInfo


def main(num_clusters, output_folder, ligand_resname, atom_ids, cpus, topology=None):
    extractCoords.main(lig_resname=ligand_resname, non_Repeat=True, atom_Ids=atom_ids, nProcessors=cpus, parallelize=False, topology=topology)
    trajectoryFolder = "allTrajs"
    trajectoryBasename = "traj*"
    stride = 1
    clusterCountsThreshold = 0
    folders = utilities.get_epoch_folders(".")
    folders.sort(key=int)

    clusteringObject = cluster.Cluster(num_clusters, trajectoryFolder,
                                       trajectoryBasename, alwaysCluster=False,
                                       stride=stride)
    clusteringObject.clusterTrajectories()
    clusteringObject.eliminateLowPopulatedClusters(clusterCountsThreshold)
    clusterCenters = clusteringObject.clusterCenters
    centersInfo = get_centers_info(trajectoryFolder, trajectoryBasename, num_clusters, clusterCenters)
    COMArray = [centersInfo[i]['center'] for i in xrange(num_clusters)]
    if output_folder is not None:
        outputFolder = os.path.join(output_folder, "")
        if not os.path.exists(outputFolder):
            os.makedirs(outputFolder)
    else:
        outputFolder = ""
    writePDB(COMArray, outputFolder+"clusters_%d_KMeans_allSnapshots.pdb" % num_clusters)
    writeInitialStructures(centersInfo, outputFolder+"initial_%d.pdb", topology=topology)

if __name__ == "__main__":
    n_clusters, lig_name, atom_id, output = parseArgs()
    main(n_clusters, output, lig_name, atom_id)
