import numpy as np
import prody as pd
import msm_pele.Box.box_helpers as mbp


def build_box(cluster_centers, env):

    # Calculate binding site point
    BS_location = mbp.get_binding_site_position(cluster_centers, env)

    # Calculate centroid
    centroid = mbp.find_centroid(cluster_centers)

    # Get center and radius along the exit path
    center, radius = get_centers(env.adap_ex_input, BS_location, centroid)

    return center, radius, "multiple"

def get_centers(system, BS_location, centroid, radius=9):
    # Equation to get the radius of one sphere if we want
    # to fit n spheres in a segment of distance d
    # radius=d/(2n-2)
    current_point = BS_location
    exit_func = lambda x: BS_location + x * (centroid-BS_location)/np.linalg.norm(centroid-BS_location)
    radius = radius if 100 > radius > 5 else 12
    alpha = 0
    # Check the distance to the centroid because we want to pass the centroid
    # and go two times that distance.
    centers = []
    while find_contacts(system, current_point):
        current_point = exit_func(alpha)
        centers.append(current_point)
        alpha += radius
        print("Point {}".format(current_point))
        print("Next distance {}".format(np.linalg.norm(BS_location-current_point)))
    radius = [radius]*len(centers)
    return centers, radius

def find_contacts(system, point):
    atoms = pd.parsePDB(system)
    contacts = pd.measure.Contacts(atoms)
    return contacts.select(5, point) 

