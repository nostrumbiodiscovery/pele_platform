#import constants

template = """
{
"waterSites":

    [
        {
        "watersToPerturb": {"links": {"ids": ["W:4"] }},
        "Box": {"radius": 6, "fixedCenter": {}, "type": "sphericalBox"}
        },
        {
        "watersToPerturb": {"links": {"ids": ["W:5", "W:6"] }},
        "Box": {"radius": 6, "fixedCenter": {}, "type": "sphericalBox"}
        }
    ]
}
"""

class PPIParams(object):


    def __init__(self, args): #, water_indexes):
        self.ppi = args.ppi
        self.center_of_interface = args.center_of_interface
        self.n_components = args.n_components if args.n_components else self.simulation_params.get("n_components", 10)
        #1) center of mass of each ligand on the input, separate 2wat, 2wat, 2wat
        #2) Buid waterstring (above)
        # template.format(com, com)
        #3) self.waters = WAT.format()
