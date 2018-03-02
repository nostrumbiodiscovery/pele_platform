import os
import sys
from MSM_PELE.Helpers import helpers, template_builder
import AdaptivePELE.adaptiveSampling as ad


class SimulationBuilder(template_builder.TemplateBuilder):

    def __init__(self, file, keywords, *args, **kwargs):

        self.file = file
        self.keywords = keywords

        self.ad_opt_new = ["false" if opt is False else opt for opt in locals()["args"]]
        self.ad_opt = ["true" if opt is True else opt for opt in self.ad_opt_new]


        self.replace = {keyword : value for keyword, value in zip(self.keywords, self.ad_opt)}

        super(SimulationBuilder, self).__init__(self.file, self.replace)

    def run(self, hook=False):
        with helpers.cd(os.path.dirname(self.file)):
            if hook:
				ad.main(self.file, clusteringHook=self.interactive_clustering)
            else:
				ad.main(self.file)

    def interactive_clustering(self, clusters, paths):
		clusters_info = clusters.__getstate__() 
		if len(clusters_info["clusters"]) == 1:
			for cluster in clusters_info["clusters"]:
				clusters_info["thresholdCalculator"].values[0] = 1.25
			return clusters, True
		elif len(clusters_info["clusters"]) >= 3 and clusters_info["thresholdCalculator"].values[0] == 1:
			clusters_info["thresholdCalculator"].values[0] = 1.75
			return clusters, True
		else:
			return clusters, False
