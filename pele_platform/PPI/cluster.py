import os
import pandas as pd
from sklearn.mixture import GaussianMixture
from pele_platform.Utilities.Helpers import bestStructs as bs
from pele_platform.Analysis.plots import _extract_coords
from multiprocessing import Pool

def cluster_best_structures(be_column, residue="LIG", topology=None, cpus=20, n_components=10, n_structs=1000):


    files_out, _, _, _, output_energy = bs.main(be_column, n_structs=n_structs, path = ".", topology=topology)
    files = []

    # find all epoch...pdb files
    for filename in os.listdir("."):
        if filename.endswith(".pdb") and filename.startswith("epoch"):
            files.append(filename)
        else:
            continue
    n_files = len(files_out)

    print("Extracting data from {} files.".format(n_files))

    snapshot = 0
    pool = Pool(cpus)
    input_pool = [[f, snapshot, residue, topology] for f in files_out]  
    all_coords = pool.map(_extract_coords, input_pool)

    # cluster
    print("Creating clusters.")
    cluster = GaussianMixture(n_components, covariance_type='full', random_state=42)
    labels = cluster.fit_predict(all_coords)

    # get lowest energy representative
    clustered_lig = pd.DataFrame(list(zip(files_out, labels, output_energy)), columns=["file_name", "cluster_ID", "binding_energy"])
    clustered_lig = clustered_lig.sort_values("binding_energy", ascending=True).groupby("cluster_ID").first()
    clustered_lig.to_csv("clustering_output.csv")

    output_files = clustered_lig["file_name"].values.tolist()
    directory = "refinement_input"

    # copy selected files to create input for refinement
    if not os.path.isdir(directory):
        os.mkdir(directory)
        for file in output_files:
            os.system("cp {} {}/.".format(file, directory))
    else:
        print(directory, "already exists!")

    return clustered_lig['file_name']
