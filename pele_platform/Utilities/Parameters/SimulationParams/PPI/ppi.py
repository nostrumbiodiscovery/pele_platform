import pele_platform.constants.constants as cs
from pele_platform.PPI.preparation import ligand_com
import pele_platform.Utilities.Helpers.helpers as hp

template = '{{"watersToPerturb": {{"links": {{"ids": "{index}" }}}}, "Box": {{"radius": 6, "fixedCenter": [{com}], "type": "sphericalBox"}}}}'

class PPIParams(object):

    def __init__(self, args, water_indices):
        self.ppi = args.ppi
        self.center_of_interface = args.center_of_interface
        self.n_components = args.n_components if args.n_components else self.simulation_params.get("n_components", 25)
        
        if args.waters and args.ppi:
       
            #For each input get me the COM of the ligand and indexes of the water inside that pdb
            coms = {}
            for inp in self.input:
                com = ligand_com(inp, args.residue)
                com_format = ['{:.10f}'.format(elem) for elem in com[0]]
                com = ", ".join(com_format)
                waters = hp.retrieve_all_waters(inp)
                coms[com] = waters
            
            water_string = []
            
            for com in coms.keys():
                for water in coms[com]:
                    water_string.append(template.format(index=water, com=com))

            self.waters = [w.strip("'") for w in water_string]
            self.water = cs.WATER_PPI.format(self.waters, self.water_temp, self.water_trials, self.water_overlap, self.water_constr).replace("'","")
            print("self.water string from PPI params", self.water)
