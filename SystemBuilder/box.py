import sys
sys.path.append("/home/dsoler/PelePlop")
import math
from Adaptive import template_builder as tb

KEYWORDS = ["CENTER_X", "CENTER_Y", "CENTER_Z", "V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8"]
COORD = "{:>11.3f}{:>8.3f}{:>8.3f}"
CENTER = "{:.3f}"


class BoxBuilder(tb.TemplateBuilder):
    """
       Create a box.pdb file with a cuadratic box
       based on a sphere grid.
    """

    def __init__(self, cm_receptor, cm_ligand, file):
        self.cm_receptor = cm_receptor
        self.cm_ligand = cm_ligand        
        self.file = file
        

    def create(self):
        """
          Description: Create pele box by
          defining center and size.
        """
        self.center = [(cm_r + cm_l)/2 for cm_r, cm_l in zip(self.cm_receptor, self.cm_ligand)]
        cm_rx, cm_ry, cm_rz = self.cm_receptor
        cm_lx, cm_ly, cm_lz = self.cm_ligand
        self.r = (math.sqrt((cm_rx-cm_lx)**2 + (cm_ry-cm_ly)**2 + (cm_rz-cm_lz)**2)/2.0)+5
        
        
    def build_box_pdb(self):
        self.build_box_crd()
        self.replace =  {keyword : value for keyword, value in zip(KEYWORDS, self.values)}
        super(BoxBuilder, self).__init__(self.file, self.replace)	

    def build_box_crd(self):
        cx, cy, cz = self.center
        self.v1 = COORD.format(cx-self.r,cy-self.r, cz-self.r)
        self.v2 = COORD.format(cx+self.r,cy-self.r, cz-self.r)
        self.v3 = COORD.format(cx+self.r,cy-self.r, cz+self.r)
        self.v4 = COORD.format(cx-self.r,cy-self.r, cz+self.r)
        self.v5 = COORD.format(cx-self.r,cy+self.r, cz-self.r)
        self.v6 = COORD.format(cx+self.r,cy+self.r, cz-self.r)
        self.v7 = COORD.format(cx+self.r,cy+self.r, cz+self.r)
        self.v8 = COORD.format(cx-self.r,cy+self.r, cz+self.r)
        self.cx = CENTER.format(cx)
        self.cy = CENTER.format(cy)
        self.cz = CENTER.format(cz)
    
        self.values = [cx, cy, cz, self.v1, self.v2, self.v3, self.v4,
                       self.v5, self.v6, self.v7, self.v8]

if __name__ == "__main__":
     BoxBuilder([11.910, 59.189, 55.414], 15, "/home/dsoler/box.pdb")
