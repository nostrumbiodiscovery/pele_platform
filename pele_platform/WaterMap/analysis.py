import glob
import numpy as np
import os
import re


class Grid:

    def __init__(self, center, radius):
        self.center = center  # user-defined for WaterMap, otherwise ligand COI
        self.radius = radius  # user-defined for WaterMap, otherwise max coordinates of the ligand
        self.side = int(radius) * 2
        self.n_voxels = int(self.side) ** 3
        self.active_voxels = []

    def generate_voxels(self):
        center = self.center
        side = self.side
        self.voxels = []

        print("Creating {} voxels.".format(self.n_voxels))

        v1 = np.add(center, [-side / 2, -side / 2, -side / 2])  # lowest coordinate vertex
        v8 = np.add(center, [side / 2, side / 2, side / 2])  # highest coordinate vertex

        first_voxel_center = v1 + np.array([0.5, 0.5, 0.5])
        last_voxel_center = v8 + np.array([-0.5, -0.5, -0.5])

        samples = self.side

        x = np.linspace(first_voxel_center[0], last_voxel_center[0], samples)
        y = np.linspace(first_voxel_center[1], last_voxel_center[1], samples)
        z = np.linspace(first_voxel_center[2], last_voxel_center[2], samples)

        x_coord, y_coord, z_coord = np.meshgrid(x, y, z)
        coordinate_grid = np.array([x_coord, y_coord, z_coord])

        voxels = list(zip(coordinate_grid[0, :, :].reshape(1, self.n_voxels).tolist()[0],
                     coordinate_grid[1, :, :].reshape(1, self.n_voxels).tolist()[0],
                     coordinate_grid[2, :, :].reshape(1, self.n_voxels).tolist()[0]))

        for v in voxels:
            voxel = Voxel(v)
            self.voxels.append(voxel)

    def add_active_voxel(self, voxel):
        self.active_voxels.append(voxel)


class Snapshot:

    def __init__(self, file, model, atom_lines):
        self.file = file
        self.model = model
        self.atom_lines = atom_lines

    def generate_properties(self):
        self.snapshot_id = "{}_{}".format(self.file, self.model)


class WaterSnapshot(Snapshot):

    def __init__(self, file, model, atom_lines):
        super().__init__(file, model, atom_lines)
        self.water_coords, self.oxygen_coords = self.extract_water_coords(self.atom_lines)
        self.orientation = self.calculate_orientation(self.water_coords)

    def extract_water_coords(self, water_lines):
        water = [line[32:56].strip().split() for line in water_lines]

        return water, water[0]  # oxygen coordinates

    def calculate_orientation(self, water_coords):

        coord_o, coord_h1, coord_h2 = water_coords
        coord_o = [float(c) for c in coord_o]
        coord_h1 = [float(c) for c in coord_h1]
        coord_h2 = [float(c) for c in coord_h2]

        v1 = np.subtract(coord_h1, coord_o)
        v2 = np.subtract(coord_h2, coord_o)

        v_middle = np.add(v1, v2)
        v_middle_magnitude = np.sqrt(v_middle[0] ** 2 + v_middle[1] ** 2 + v_middle[2] ** 2)

        orientation = [np.degrees(np.arccos(v_middle[0] / v_middle_magnitude)),
                       np.degrees(np.arccos(v_middle[1] / v_middle_magnitude)),
                       np.degrees(np.arccos(v_middle[2] / v_middle_magnitude))]

        return orientation

    def check_voxel(self, voxels):

        distances = {}

        for voxel in voxels:
            oxygen_coords = [float(c) for c in self.oxygen_coords]
            dist = np.subtract(oxygen_coords, voxel.voxel_center)
            dist = np.linalg.norm(dist)
            distances[voxel] = dist

        min_dist = min(list(distances.values()))
        result = [key for key, value in distances.items() if value == min_dist]

        self.voxel = result[0]

        return self.voxel


class Voxel:

    def __init__(self, v):
        self.voxel_center = v
        self.snapshots = []

    def add_snapshot(self, snap):
        self.snapshots.append(snap)

    def calculate_entropy(self):
        pass

    def calculate_enthalpy(self):
        pass


def extract_snapshots(output_simulation):

    print("Extracting water snapshots from {}.".format(output_simulation))
    all_files = glob.glob(os.path.join(output_simulation, '**', 'trajectory*'))
    all_snapshots = []

    if len(all_files) == 0:
        print("No output files detected.")
    else:
        for f in all_files:
            with open(f, "r") as file:
                file_content = file.read()
                matches = re.findall(r'(MODEL.+?ENDMDL)', file_content, re.DOTALL)

                for m in matches:
                    model = re.search(r'MODEL\s+(\d+)', m).group(1)
                    if model != '1':
                        water_match = re.findall(r'.+HOH.+', m)
                        snapshot = WaterSnapshot(f, model, water_match)
                        all_snapshots.append(snapshot)

    return all_snapshots


def main(center, radius, output_simulation):

    grid = Grid(center, radius)
    grid.generate_voxels()
    snapshots = extract_snapshots(output_simulation)

    for s in snapshots:
        s.generate_properties()
        voxel = s.check_voxel(grid.voxels)
        voxel.add_snapshot(s)
        grid.add_active_voxel(voxel)


if __name__ == "__main__":

    center = [29.842, 48.919, 61.586]
    radius = 20
    output_simulation = "water_sim"

    grid = Grid(center, radius)
    grid.generate_voxels()
    snapshots = extract_snapshots(output_simulation)

    for s in snapshots:
        s.generate_properties()
        voxel = s.check_voxel(grid.voxels)
        voxel.add_snapshot(s)
        grid.add_active_voxel(voxel)



