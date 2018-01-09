import os
import sys
import subprocess
from Helpers.pele_env import cd
from template_builder import TemplateBuilder
from AdaptivePELE.adaptiveSampling import main


class AdaptiveBuilder(TemplateBuilder):

    def __init__(self, file, keywords, *args, **kwargs):
        
        self.file = file
        self.keywords = keywords

        self.ad_opt_new = ["false" if opt == False else opt for opt in locals()["args"]]
        self.ad_opt = ["true" if opt == True else opt for opt in self.ad_opt_new]

        
        self.replace =  {keyword : value for keyword, value in zip(self.keywords, self.ad_opt)}

        super(AdaptiveBuilder, self).__init__(self.file, self.replace)


    def run(self):

    	with cd(os.path.dirname(self.file)):
    	 main(self.file)








