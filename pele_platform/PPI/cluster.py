import os
import pandas as pd
from sklearn.mixture import GaussianMixture
from pele_platform.Utilities.Helpers import bestStructs, helpers


def cluster_best_structures(
    be_column: int,
    residue="LIG",
    topology=None,
    cpus=20,
    n_components=10,
    n_structs=1000,
    directory=".",
    logger=None,
):

    files_out, _, _, _, output_energy = bestStructs.main(
        be_column, n_structs=n_structs, path=".", topology=topology, logger=logger
    )
    n_files = len(files_out)
    logger.info("Extracting data from {} files.".format(n_files))

    snapshot = 0
    input_pool = [[f, snapshot, residue, topology] for f in files_out]
    all_coords = helpers.parallelize(_extract_coords, input_pool, 1)

    # cluster
    logger.info("Creating clusters.")
    cluster = GaussianMixture(n_components, covariance_type="full", random_state=42)
    labels = cluster.fit_predict(all_coords)

    # get lowest energy representative
    clustered_lig = pd.DataFrame(
        list(zip(files_out, labels, output_energy)),
        columns=["file_name", "cluster_ID", "binding_energy"],
    )
    clustered_lig = (
        clustered_lig.sort_values("binding_energy", ascending=True)
        .groupby("cluster_ID")
        .first()
    )
    clustered_lig.to_csv("clustering_output.csv")

    output_files = clustered_lig["file_name"].values.tolist()
    directory = os.path.join(directory, "refinement_input")

    # copy selected files to create input for refinement
    if not os.path.isdir(directory):
        os.makedirs(directory, exist_ok=True)
    for file in output_files:
        os.system("cp {} {}/.".format(file, directory))

    return clustered_lig["file_name"]


def _extract_coords(info):
    import numpy as np
    import mdtraj

    p, v, resname, topology = info
    # Most time consuming step 0.1
    traj = mdtraj.load_frame(p, v, top=topology)
    atoms_info = traj.topology.to_dataframe()[0]
    condition = atoms_info["resName"] == resname
    atom_numbers_ligand = atoms_info[condition].index.tolist()
    coords = []
    for atom_num in atom_numbers_ligand:
        try:
            coords.extend(traj.xyz[0][atom_num].tolist())
        except IndexError:
            continue
    return np.array(coords).ravel() * 10
