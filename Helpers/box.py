# from matplotlib import pyplot
# from mpl_toolkits.mplot3d import Axes3D
import os
import random
from scipy.spatial import distance
import numpy as np
import argparse
import math
from operator import itemgetter, attrgetter
import MSM_PELE.Helpers.best_structs as best_structs
import MSM_PELE.Helpers.template_builder as tb

__author__ = "Daniel Soler Viladrich"
__email__ = "daniel.soler@nostrumbiodiscovery.com"

KEYWORDS = ["CENTER_X", "CENTER_Y", "CENTER_Z", "V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8"]
COORD = "{:>11.3f}{:>8.3f}{:>8.3f}"
CENTER = "{:.3f}"


def parseargs():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('bs', type=int, nargs='+', help='an integer for the accumulator')
    parser.add_argument('--points', '-p', type=str, help='File where to find all coord')
    args = parser.parse_args()
    return args.bs, args.points


def main(path, clusters, bs):
    # fig = pyplot.figure()
    # ax = Axes3D(fig)

    sasa_info = best_structs.main(path)    
    sasa_points = get_sasa_points(path, sasa_info)
    points = get_points(clusters)
    centroid = find_centroid(points)

    angle_points = find_angle_lenght(bs, centroid, sasa_points)
    min_angl_points = sorted(angle_points, key=itemgetter(1))
    chosen_point = min_angl_points[0][0]

    radius = (distance.euclidean(bs, chosen_point) / 2.0) + 4
    center = [(final + initial) / 2.0 for initial, final in zip(bs, chosen_point)]
  #  ax.scatter(*decompose(points))
  #  ax.plot_wireframe(*WireframeSphere(center, radius), color="b", alpha=0.5)
  #  ax.scatter(*bs, c='red')
  #  ax.scatter(*chosen_point, c='red')
 #   ax.scatter(*centroid, c='red')

    #pyplot.show()

    return center, radius

def get_sasa_points(path, sasa_info):
    points = []
    for _, info in sasa_info.items():
        epoch, report, value, model = info
        coord_file = os.path.join(path, "{}/extractedCoordinates/coord_{}.dat".format(epoch, report))
        with open(coord_file, 'r') as f:
                lines = f.readlines()
        coord = [ float(crd) for crd in lines[model].split()[1:]]
        points.append(coord)
    return points




def get_points(pdb):
    with open(pdb, 'r') as f:
        lines = [line.split() for line in f if line.startswith("HETATM")]
        points = [[float(line[6]), float(line[7]), float(line[8])] for line in lines]
        return points




def find_centroid(points):
    x, y, z = decompose(points)
    n_points = len(points)
    centroid = (sum(x) / n_points, sum(y) / n_points, sum(z) / n_points)
    return centroid


def decompose(points):
    crd_x = [ x for x, y, z in points]
    crd_y = [ y for x, y, z in points]
    crd_z = [ z for x, y, z in points]
    return crd_x, crd_y, crd_z


def find_angle_lenght(bs, centroid, points):
    point_angle_lenght = []
    bs_centroid = [final - initial for initial, final in zip(bs, centroid)]
    point_centroid = []
    for point in points:
        point_centroid.append([final - initial for initial, final in zip(centroid, point)])
    for point, vector in zip(points, point_centroid):
        scalar = 0
        for initial, final in zip(vector, bs_centroid):
            scalar += initial * final
        magnitude_bs_cent = math.sqrt((bs_centroid[0])**2 + (bs_centroid[1])**2 + (bs_centroid[2])**2)
        magnitude_point_cent = math.sqrt((vector[0])**2 + (vector[1])**2 + (vector[2])**2)
        angle = math.acos(scalar / (magnitude_point_cent * magnitude_bs_cent))
        point_angle_lenght.append([point, angle, distance.euclidean(centroid, point)])
    return point_angle_lenght



def WireframeSphere(centre, radius, n_meridians=20, n_circles_latitude=None):
    """
    Create the arrays of values to plot the wireframe of a sphere.

    """
    if n_circles_latitude is None:
        n_circles_latitude = max(n_meridians / 2, 4)
    u, v = np.mgrid[0:2 * np.pi:n_meridians * 1j, 0:np.pi:n_circles_latitude * 1j]
    sphere_x = centre[0] + radius * np.cos(u) * np.sin(v)
    sphere_y = centre[1] + radius * np.sin(u) * np.sin(v)
    sphere_z = centre[2] + radius * np.cos(v)
    return sphere_x, sphere_y, sphere_z

def build_box(center, radius, file):

    cx, cy, cz = center
    v1 = COORD.format(cx - radius, cy - radius, cz - radius)
    v2 = COORD.format(cx + radius, cy - radius, cz - radius)
    v3 = COORD.format(cx + radius, cy - radius, cz + radius)
    v4 = COORD.format(cx - radius, cy - radius, cz + radius)
    v5 = COORD.format(cx - radius, cy + radius, cz - radius)
    v6 = COORD.format(cx + radius, cy + radius, cz - radius)
    v7 = COORD.format(cx + radius, cy + radius, cz + radius)
    v8 = COORD.format(cx - radius, cy + radius, cz + radius)
    cx = CENTER.format(cx)
    cy = CENTER.format(cy)
    cz = CENTER.format(cz)

    values = [cx, cy, cz, v1, v2, v3, v4, v5, v6, v7, v8]

    replace = {keyword: value for keyword, value in zip(KEYWORDS, values)}

    tb.TemplateBuilder(file, replace)

    return file

if __name__ == '__main__':
        #with pele.cd("/home/dsoler/1sqa/output_adaptive_short/"):
        #    cl.main(num_clusters=40, output_folder="/home/dsoler/1sqa/output_clustering/", ligand_resname="UI1", atom_ids="")
       	center, radius = main("/scratch/jobs/dsoler/test/STR_Pele/output_adaptive_exit" ,"/scratch/jobs/dsoler/test/STR_Pele/output_clustering/clusters_40_KMeans_allSnapshots.pdb" , [-20.332, 59.897, 2.8323])
        box = build_box(center, radius, "/scratch/jobs/dsoler/test/STR_Pele/box.pdb")

