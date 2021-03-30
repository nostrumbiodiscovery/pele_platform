from Bio.PDB import PDBParser
import glob, numpy, os
import pandas as pd

parser = PDBParser()


def get_atom_from_ref(file, atomname='C29', resname="FMM"):
    reference = parser.get_structure("reference", file)
    for residue in reference.get_residues():
        if residue.resname == resname:
            for atom in residue.get_atoms():
                if atom.name == atomname:
                    atom_ref = atom
                    break
            break
    return atom_ref


def get_best_distance(pdb_file, reference_point, resname="GRW"):
    """
    Finds fragment atom closest to the user-defined reference point.

    Parameters
    ----------
    pdb_file : str
        Path to PDB file.
    reference_point : list[float]
        Coordinates of the reference point to which the distance should be calculated.
    resname : str
        Residue name of the grown fragment, default = "GRW".
    Returns
    -------
        Distance from the closest atom to the reference point and closes atom name.
    """
    struct = parser.get_structure("epoch", pdb_file)
    ref_vector = numpy.array(reference_point)

    best_dist = None

    for residue in struct.get_residues():
        if residue.resname == resname:
            for atom in residue.get_atoms():
                atom_vector = numpy.array(atom.get_coord())
                dist = numpy.linalg.norm(atom_vector - ref_vector)

                if not best_dist or dist < best_dist:
                    best_atom = atom
                    best_dist = dist

    return best_dist, best_atom


def get_distance(reference_point, path, resname):
    best_dist = None

    for file in glob.glob(path):
        distance, atom = get_best_distance(file, reference_point, resname=resname)
        if not best_dist or distance < best_dist:
            best_dist = distance
            best_file_distance = file
            best_atom = atom

    return best_dist, best_file_distance, best_atom


def get_distance_be(reference_point, path, resname):
    best_dist = ''
    best_file_distance = ''
    best_atom = ''
    best_be = 0
    best_file_be = ''
    dist_list = []
    be_list = []

    for file in glob.glob(path):
        distance, atom = get_best_distance(file, reference_point, resname=resname)
        dist_list.append(distance)
        if not best_dist or distance < best_dist:
            best_dist = distance
            best_file_distance = file
            best_atom = atom
        be = file.split("_BindingEnergy")[1]
        be = float(be[:-4])
        be_list.append(be)
        if be == 0 or be < best_be:
            best_be = be
            best_file_be = file

    dist_list_norm, be_list_norm = normalize_lists(dist_list, be_list)
    norm_list = numpy.sqrt(numpy.multiply(numpy.power(dist_list_norm, 2), numpy.power(be_list_norm, 2)))
    best_normalize = min(norm_list)
    index_best_norm = numpy.where(norm_list == best_normalize)
    best_norm_file = glob.glob(path)[index_best_norm[0][0]]
    dist_best_norm = dist_list[index_best_norm[0][0]]
    be_best_norm = be_list[index_best_norm[0][0]]
    return best_dist, best_file_distance, best_atom, best_be, best_file_be, best_normalize, best_norm_file, dist_best_norm, be_best_norm


def normalize_lists(distance_list, be_list):
    max_dist = max(distance_list)
    distance_list = [x / max_dist for x in distance_list]
    min_be = min(be_list)
    be_list = numpy.power([x / min_be for x in be_list], -1)
    return distance_list, be_list


def main(path, atom_coords, resname="GRW", pattern='', reference_file=""):
    if not atom_coords:
        atom_coords = get_atom_from_ref(reference_file).get_coord()

    search_path = os.path.join(path, "*{}*".format(os.path.splitext(pattern)[0]))
    config_files = glob.glob(search_path)

    files = []
    best_files_found_distance = []
    dists = []
    best_files_found_be = []
    bes = []
    norm = []
    best_files_found_norm = []
    dist_norm_list = []
    be_norm_list = []

    for file in config_files:
        if os.path.exists(os.path.join(file, "top_result")):
            path = os.path.join(file, "top_result", "epoch*")
            best_dist, best_file_distance, best_atom, best_be, best_file_be, best_normalize, best_norm_file, dist_best_norm, be_best_norm = get_distance_be(
                atom_coords, path, resname)
            dist_math = format(best_dist, '.15g')

            files.append(file)
            best_files_found_distance.append(best_file_distance)
            dists.append(dist_math)
            best_files_found_be.append(best_file_be)
            bes.append(best_be)
            norm.append(best_normalize)
            best_files_found_norm.append(best_norm_file)
            dist_norm_list.append(dist_best_norm)
            be_norm_list.append(be_best_norm)

    dataframe = {'File': files, 'BestFileDistance': best_files_found_distance, 'Distance': dists,
                 'BestFileBE': best_files_found_be, 'BE': bes, 'BestFileNormalization': best_files_found_norm,
                 'BestNormalization': norm, 'DistanceBestNormalization': dist_norm_list,
                 'BEBestNormalization': be_norm_list}
    df = pd.DataFrame(dataframe,
                      columns=['File', 'BestFileDistance', 'Distance', 'BestFileBE', 'BE', 'BestFileNormalization',
                               'BestNormalization', 'DistanceBestNormalization', 'BEBestNormalization'])
    df.to_csv('point_analysis.csv', index=False, header=True)
