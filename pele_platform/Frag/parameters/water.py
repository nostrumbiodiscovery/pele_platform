from pele_platform.Utilities.Helpers import water as wt


class FragWaterParams():

    def __init__(self, pdb):
        self.pdb = pdb
        self.initialize_waters()

    def initialize_waters(self):
        self.water_object = wt.WaterIncluder([self.pdb], self.n_waters,
            user_waters=self.waters, ligand_perturbation_params=self.parameters,
            water_center=self.water_center, water_radius=self.water_radius,
            allow_empty_selectors=self.allow_empty_selectors, water_temp=self.water_temp,
            water_trials=self.water_trials, water_overlap=self.water_overlap,
            water_constr=self.water_constr, test=self.test, water_freq=self.water_freq,
            ligand_residue=self.residue)
        self.water_object.run()
    
