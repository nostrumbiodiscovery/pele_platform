import os
import sys
from template_builder import TemplateBuilder


class AdaptiveBuilder(TemplateBuilder):

    def __init__(self, file, keywords, *args, **kwargs):
        
        self.file = file
        self.keywords = keywords

        self.ad_opt = [None if opt == 'None' else opt for opt in locals()["args"]]
        self.ad_opt_new = [False if opt == 'False' else opt for opt in self.ad_opt]
        self.ad_opt = [True if opt == 'True' else opt for opt in self.ad_opt_new]

        
        self.replace =  {keyword : value for keyword, value in zip(self.keywords, self.ad_opt)}

        super(AdaptiveBuilder, self).__init__(self.file, self.replace)


    def run(self):

    	cmd = "python -m AdaptivePELE.adaptiveSampling {}".format(self.file)

    	print(cmd)

    	#process = subprocess.Popen(cmd,shell=False,stdin=None,stdout=None,stderr=None,close_fds=True, preexec_fn=os.setsid)








