import numpy as np
import msm_pele.Analysis.plotMSMAdvancedInfo as pm


def get_points(pdb):
    with open(pdb, 'r') as f:
        lines = [line.split() for line in f if line.startswith("HETATM")]
        points = [[float(line[6]), float(line[7]), float(line[8])] for line in lines]
        return points


def find_centroid(points):
    x = [cx for cx, cy, cz in points]
    y = [cy for cx, cy, cz in points]
    z = [cz for cx, cy, cz in points]
    n_points = len(points)
    centroid = (sum(x) / n_points, sum(y) / n_points, sum(z) / n_points)
    return centroid


def get_binding_site_position(cluster_centers, env):
    minpos = pm.get_min_Pos(env.system_fix, env.residue)
    distances = np.linalg.norm(cluster_centers - minpos, axis=1)
    min_distance = np.argmin(distances)
    BS_location = cluster_centers[min_distance]
    return BS_location
