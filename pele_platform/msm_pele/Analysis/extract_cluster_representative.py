import numpy as np
import shutil
import os
import argparse
import msm_pele.Analysis.plotMSMAdvancedInfo as pt



def main(native, resname, destFolder, dg1, dg2, dist1, dist2, output="most_important_clusters"):
    minPos = pt.get_min_Pos(native, resname)
    clusters = np.loadtxt(os.path.join(destFolder, "clusterCenters_%d.dat" % i))
    distance = np.linalg.norm(clusters-minPos, axis=1)
    data = np.loadtxt(os.path.join(destFolder, "pmf_xyzg_%d.dat" % i))
    deltag = data[:, -1]
    cluster_idx = [j for j, (g, d) in enumerate(zip(deltag, distance)) if (dg2 > g > dg1) and (dist2 > d > dist1)]
    folder = os.path.join(os.path.dirname(os.path.abspath(destFolder)), "representative_structures_pdb_0_{}".format(i))
    structures_clusters = [os.path.join(folder, "cluster_{}.pdb".format(idx)) for idx in cluster_idx]
    if not os.path.exists(output):
        os.mkdir(output)
    for f in structures_clusters:
        shutil.copy(f, output)


def add_args(parser):
    parser.add_argument("native", type=str, help="Initial structure")
    parser.add_argument("resname", type=str, help="Residue name of your ligand")
    parser.add_argument("destFolder", type=str, help="MSM folder")
    parser.add_argument("--dgmin", type=float, help="DG min limit")
    parser.add_argument("--dgmax", type=float, help="DG ax limit")
    parser.add_argument("--distmin", type=float, help="dist min limit")
    parser.add_argument("--distmax", type=float, help="dist max limit")
    parser.add_argument("--output", type=str, help="output folder", default = "most_important_clusters")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Extract the representative structures of the most important clusters')
    add_args(parser)
    args = parser.parse_args()
    main(args.native, args.resname, args.destFolder, args.dgmin, args.dgmax, args.distmin, args.distmax, args.output)
