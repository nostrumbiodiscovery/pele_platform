import glob
import numpy as np
import os
import pandas as pd
import re
from sklearn.mixture import GaussianMixture

class Watermap:

    def __init__(self, output_simulation, args=None):
        self.output_simulation = output_simulation
        self.radius = 8
        self.water_center = [50, 50, 50]

    def run(self):
        self.all_snapshots = self.extract_snapshots()
        self.grid = Grid(self.water_center, self.radius)
        self.cluster_waters(self.all_snapshots)
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

        for o in orientation:
            o_coords = o[2][0]
            is_in = is_in_voxel(o_coords, self.voxels)
            o.extend(is_in)

        clusters = None
       # all_coords = [o[2] for o in orientation]
       # all_coords = [a[0]+a[1]+a[2] for a in all_coords]
       # n_components = int(len(all_coords)/10)
       # clusters = GaussianMixture(n_components, covariance_type='full', random_state=42)
       # labels = clusters.fit_predict(all_coords)
       # clusters = zip(orientation, labels) 
       # results = pd.DataFrame(list(clusters))
       # results.to_csv("clustering_output.csv")

        return clusters

    def calculate_entropy(self):
        pass

    def calculate_enthalpy(self):
        pass


class Grid:

    def __init__(self, center, radius):
        self.center = center
        self.side = radius/2
        self.vertices = self.calculate_vertices()
        self.voxels = self.calculate_voxels(self.vertices, self.side)

    def calculate_vertices(self):
        center = self.center
        side = self.side

        vertices = []
        v1 = np.add(center,[-side/2, -side/2, -side/2])
        v2 = np.add(center, [side/2, -side/2, -side/2])
        v3 = np.add(center, [-side/2, -side/2, side/2])
        v4 = np.add(center, [-side/2, side/2, -side/2])
        v5 = np.add(center, [side/2, side/2, -side/2])
        v6 = np.add(center, [side/2, -side/2, side/2])
        v7 = np.add(center, [-side/2, side/2, side/2])
        v8 = np.add(center, [side/2, side/2, side/2])
        vertices.extend([v1, v2, v3, v4, v5, v6, v7, v8])

        return vertices

    def calculate_voxels(self, vertices, side):

        self.n_voxels = int(side)**3
        print("Creating {} voxels.".format(self.n_voxels))

        first_voxel_center = vertices[0] + np.array([0.5, 0.5, 0.5])
        last_voxel_center = vertices[7] + np.array([-0.5, -0.5, -0.5])

        samples =  int(max(last_voxel_center) - min(first_voxel_center) + 1)  

        x = y = z = np.linspace(min(first_voxel_center), max(last_voxel_center), samples)
        x_coord, y_coord, z_coord = np.meshgrid(x, y, z)
        coordinate_grid = np.array([x_coord, y_coord, z_coord])
        voxels = zip(coordinate_grid[0,:,:].reshape(1, self.n_voxels).tolist()[0], coordinate_grid[1,:,:].reshape(1, self.n_voxels).tolist()[0], coordinate_grid[2,:,:].reshape(1, self.n_voxels).tolist()[0])
        self.voxels = voxels
        return voxels

def is_in_voxel(self, o_coord, voxels):
    distances = []
    for v in voxels:
        print(v)
        dist = numpy.linalg.norm(v - o_coord)
        distances.append([dist])

    in_voxel = min(distances)
    return in_voxel
                        

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


   # grid = Grid(center, radius)
   # pixels = grid.generate_pixels()
   # snapshots = find_snapshots()
   # for s in snapshots: #Parallel!!!!
   #     snapshot = Snapshot(s)
   #     snapshot.generate_properties()
   #     for pixel in pixels:
   #         if pixel.is_inside(snapshot.oxigen_coords):
   #             pixel.add(snapshot)
