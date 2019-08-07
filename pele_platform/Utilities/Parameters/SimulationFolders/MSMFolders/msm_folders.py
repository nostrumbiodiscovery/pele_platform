import os

class MSMPaths(object):


    def __init__(self, args):
        self.cluster_output = os.path.join(self.pele_dir, "output_clustering")
        self.clusters_output = os.path.join(self.cluster_output, "clusters_{}_KMeans_allSnapshots.pdb".format(self.clusters))
