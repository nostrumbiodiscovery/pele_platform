import os
import sys
import subprocess
from Helpers import helpers, template_builder
from template_builder import TemplateBuilder
import AdaptivePELE.adaptiveSampling as ad


class SimulationBuilder(template_builder.TemplateBuilder):

    def __init__(self, file, keywords, *args, **kwargs):
        
        self.file = file
        self.keywords = keywords

        self.ad_opt_new = ["false" if opt == False else opt for opt in locals()["args"]]
        self.ad_opt = ["true" if opt == True else opt for opt in self.ad_opt_new]

        
        self.replace =  {keyword : value for keyword, value in zip(self.keywords, self.ad_opt)}

        super(SimulationBuilder, self).__init__(self.file, self.replace)


    def run(self):

    	with helpers.cd(os.path.dirname(self.file)):
    	 ad.main(self.file)








