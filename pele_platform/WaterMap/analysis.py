import glob
import numpy as np
import os
import pandas as pd
import re
from sklearn.mixture import GaussianMixture


class Watermap:

    def __init__(self, output_simulation):
        self.output_simulation = output_simulation

    def run(self):
        all_snapshots = self.extract_snapshots()
        clusters = self.cluster_waters(all_snapshots)
        self.calculate_entropy()
        self.calculate_enthalpy()

    def extract_snapshots(self):
        
        print("Extracting water snapshots from {}.".format(self.output_simulation))
        all_files = glob.glob(os.path.join(self.output_simulation, '**', 'trajectory*'))
        all_snapshots = []

        if len(all_files) == 0:
            print("No output files detected.")
        else:
            for f in all_files:
                with open(f, "r") as file:
                    file_content = file.read()
                    matches = re.findall('(MODEL.+?ENDMDL)', file_content, re.DOTALL)

                    for m in matches:
                        model = re.search('MODEL\s+(\d+)', m).group(1)
                        if model != '1':
                            water_match = re.findall('.+HOH.+', m)
                            all_snapshots.append([f, model, water_match])
        return all_snapshots

    def cluster_waters(self, snapshots):
        
        water_coords = extract_water_coords(snapshots)
        orientation = calculate_orientation(water_coords)
        all_coords = [o[2] for o in orientation]
        all_coords = [a[0]+a[1]+a[2] for a in all_coords]
        n_components = int(len(all_coords)/10)
        clusters = GaussianMixture(n_components, covariance_type='full', random_state=42)
        labels = clusters.fit_predict(all_coords)
        clusters = zip(orientation, labels) 
        
        results = pd.DataFrame(list(clusters))
        results.to_csv("clustering_output.csv")

        return clusters

    def calculate_entropy(self):
        pass

    def calculate_enthalpy(self):
        pass


def extract_water_coords(snapshots):
        
    water_coords = []

    for snap in snapshots:
        file, model, water_lines = snap
        file_id = "traj{}_model{}".format(os.path.basename(file).replace(".pdb","").replace("trajectory_",""), model)
        water = [line[32:56].strip().split() for line in water_lines]
        water_coords.append([file_id, model, water])
    return water_coords


def calculate_orientation(water_coords):

    output = []

    for water in water_coords:
        file_id, model, coords = water    
        coord_o, coord_h1, coord_h2 = coords
        coord_o = [float(c) for c in coord_o]
        coord_h1 = [float(c) for c in coord_h1]
        coord_h2 = [float(c) for c in coord_h2]

        v1 = np.subtract(coord_h1,coord_o)
        v2 = np.subtract(coord_h2,coord_o)
        v_middle = np.add(v1, v2)
        v_middle_magnitutde = np.sqrt(v_middle[0]**2+v_middle[1]**2+v_middle[2]**2)
        orientation = [np.degrees(np.arccos(v_middle[0]/v_middle_magnitutde)), np.degrees(np.arccos(v_middle[1]/v_middle_magnitutde)), np.degrees(np.arccos(v_middle[2]/v_middle_magnitutde))]
        water.extend([orientation])
        output.append(water)
    return output


if __name__ == "__main__":
    analysis = Watermap("water_sim/output/")
    analysis.run()
