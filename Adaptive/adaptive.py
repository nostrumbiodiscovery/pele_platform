import os
import sys
from template_builder import TemplateBuilder


class AdaptiveBuilder(TemplateBuilder):

    def __init__(self, file, keywords, *args, **kwargs):
        
        self.file = file
        self.keywords = keywords

        ad_opt = [None if opt == 'None' else opt for opt in locals()["args"]]
        ad_opt_new = [False if opt == 'False' else opt for opt in ad_opt]
        ad_opt = [True if opt == 'True' else opt for opt in ad_opt_new]

        print(ad_opt)
        
        self.replace =  {keyword : value for keyword, value in zip(self.keywords, ad_opt)}

        super(AdaptiveBuilder, self).__init__(self.file, self.replace)




