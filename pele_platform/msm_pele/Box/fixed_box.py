import numpy as np
import prody as pd
import msm_pele.Helpers.center_of_mass as cm
import msm_pele.Box.box_helpers as mbp


def build_box(system, clusters, ligand):
    bs = cm.center_of_mass(ligand)
    points = mbp.get_points(clusters)
    centroid = mbp.find_centroid(points)
    center = find_non_contact_points(system, centroid, bs)
    radius = (np.linalg.norm(np.array(bs) - np.array(center)) + 10)
    return [center, ], [radius, ], "fixed"


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

    point = np.array(bs, dtype=float)
    number_of_contacts = True
    while number_of_contacts:
        number_of_contacts = contacts.select(5, point)
        point = np.array(point, dtype=float) + directior_unitary
    return point.tolist()
