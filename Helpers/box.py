import os
import prody as pd
from scipy.spatial import distance
import shutil
import numpy as np
import argparse
import math
import MSM_PELE.Helpers.best_structs as best_structs
import MSM_PELE.Helpers.template_builder as tb
import MSM_PELE.Helpers.helpers as hp
import MSM_PELE.Helpers.center_of_mass as cm

__author__ = "Daniel Soler Viladrich"
__email__ = "daniel.soler@nostrumbiodiscovery.com"

# BOX CONSTANTS
KEYWORDS = ["RADIUS", "CENTER_X", "CENTER_Y", "CENTER_Z", "V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8"]
COORD = "{:>11.3f}{:>8.3f}{:>8.3f}"
CENTER = "{:.3f}"


def parseargs():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('bs', type=int, nargs='+', help='an integer for the accumulator')
    parser.add_argument('--points', '-p', type=str, help='File where to find all coord')
    args = parser.parse_args()
    return args.bs, args.points


def create_box(args, env):
    BS_sasa_min, BS_sasa_max = is_exit_finish(env.adap_ex_output, args.test)
    if args.box:
        center, radius = retrieve_box_info(args.box, env.clusters_output)
        shutil.copy(args.box, os.path.join(env.pele_dir, "box.pdb"))
    else:
        center_mass = cm.center_of_mass(env.ligand_ref)
        center, radius = build_box(env.adap_ex_input, env.clusters_output, center_mass)
        if args.user_center and args.user_radius:
            center = args.user_center
            radius = args.user_radius
        box_to_pdb(center, radius, env.box_temp)
    return center, radius, BS_sasa_min, BS_sasa_max


def build_box(system, clusters, bs):

    points = get_points(clusters)
    centroid = find_centroid(points)
    center = find_non_contact_points(system, centroid, bs) 

    radius = (distance.euclidean(bs, center) + 10) 

    remove_clusters_out_of_box(os.path.dirname(clusters), center, radius, points)

    return center, radius

def remove_clusters_out_of_box(cluster_directory, center, radius, points):
    cx, cy, cz = center
    for i, point in enumerate(points):
        x,y,z = point
        if ((x-cx)**2 + (y-cy)**2 + (z-cz)**2) > (radius**2):
            try:
                os.remove(os.path.join(cluster_directory, "initial_{}.pdb".format(i)))
            except OSError:
                pass


def find_non_contact_points(system, centroid, bs):
    """
        Find 0 contact point with protein
        in the direction produced by the 
        points of the binding side and the 
        centroide.
    """

    direction = np.array(centroid, dtype=float) - np.array(bs, dtype=float)
    directior_unitary = direction / np.linalg.norm(direction)
    atoms = pd.parsePDB(system)
    contacts = pd.measure.Contacts(atoms)

    point=np.array(bs, dtype=float)
    number_of_contacts = False
    while number_of_contacts > 0 or number_of_contacts is False:
        number_of_contacts = contacts. select(5, point)
        point = np.array(point, dtype=float) + directior_unitary
        print(point)
        print(number_of_contacts)
    return point.tolist()

def get_sasa_points(path, max_sasa_structs):
    points = []
    for _, info in max_sasa_structs.items():
        epoch, report, value, model = info
        coord_file = os.path.join(path, "{}/extractedCoordinates/coord_{}.dat".format(epoch, report))
        with open(coord_file, 'r') as f:
                lines = f.readlines()
        coord = [float(crd) for crd in lines[model].split()[1:]]
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
    crd_x = [x for x, y, z in points]
    crd_y = [y for x, y, z in points]
    crd_z = [z for x, y, z in points]
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

def retrieve_box_info(box, clusters):
    with open(box, 'r') as f:
        lines = hp.preproces_lines(f.readlines()) 
        try:
            center = [float(coord) for coord in lines[1][5:9]]
        except ValueError:
            raise ValueError("{} not valid. Check the file is not a template.")
        radius = float(lines[2][2])
    print(center, radius)
    points = get_points(clusters)
    remove_clusters_out_of_box(os.path.dirname(clusters), center, radius, points)
    return center, radius

def box_to_pdb(center, radius, file):

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

    values = [radius, cx, cy, cz, v1, v2, v3, v4, v5, v6, v7, v8]

    replace = {keyword: value for keyword, value in zip(KEYWORDS, values)}

    tb.TemplateBuilder(file, replace)


def is_exit_finish(path, test):
    return best_structs.main(path, test=test)
