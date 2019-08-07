import os

class GlidePaths(object):


    def __init__(self, args):
        self.glide_template = os.path.join(self.pele_dir, "glide.in") 
        self.glide_structs = os.path.join(self.pele_dir, "glide_calculations", "structures")
