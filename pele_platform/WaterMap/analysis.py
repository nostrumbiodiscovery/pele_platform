import glob
import numpy as np
import os
import re


class Grid:

    def __init__(self, box_center, box_radius):
        """
        Initializes the Grid class.

        Parameters
        ----------
        box_center : List[str]
            Box center defined by the user, otherwise ligand center of mass.
        box_radius : float
            Box radius defined by the user, otherwise maximum coordinates of the ligand.
        """
        self.center = [float(number) for number in box_center]
        self.radius = box_radius
        self.side = int(box_radius) * 2
        self.n_voxels = int(self.side) ** 3
        self.active_voxels = []
        self.voxels = list()

    def generate_voxels(self):
        """
        Generates a grid of 1A around.
        """
        # lowest coordinate vertex
        v1 = np.add(
            self.center, [-self.radius, -self.radius, -self.radius]
        )

        # highest coordinate vertex
        v8 = np.add(
            self.center, [self.radius, self.radius, self.radius]
        )

        first_voxel_center = v1 + np.array([0.5, 0.5, 0.5])
        last_voxel_center = v8 + np.array([-0.5, -0.5, -0.5])

        samples = self.side  # number of voxels per side

        x = np.linspace(first_voxel_center[0], last_voxel_center[0], samples)
        y = np.linspace(first_voxel_center[1], last_voxel_center[1], samples)
        z = np.linspace(first_voxel_center[2], last_voxel_center[2], samples)

        x_coord, y_coord, z_coord = np.meshgrid(x, y, z)
        coordinate_grid = np.array([x_coord, y_coord, z_coord])

        voxels = list(
            zip(
                coordinate_grid[0, :, :].reshape(1, self.n_voxels).tolist()[0],
                coordinate_grid[1, :, :].reshape(1, self.n_voxels).tolist()[0],
                coordinate_grid[2, :, :].reshape(1, self.n_voxels).tolist()[0],
            )
        )

        for v in voxels:
            voxel_object = Voxel(v)
            self.voxels.append(voxel_object)

    def add_active_voxel(self, voxel):
        self.active_voxels.append(voxel)


class Snapshot:
    def __init__(self, file, model, atom_lines):
        """
        Initializes Snapshot class.

        Parameters
        ----------
        file : str
            Path to file (trajectory).
        model : int
            Number of model in a trajectory.
        atom_lines : List[str]
            PDB lines referring to the atom.
        """
        self.file = file
        self.model = model
        self.atom_lines = atom_lines
        self.snapshot_id = None

    def generate_properties(self):
        """
        Generates properties:
            - snapshot ID from file name and model number.
        """
        self.snapshot_id = "{}_{}".format(self.file, self.model)


class WaterSnapshot(Snapshot):

    def __init__(self, file, model, atom_lines):
        """
        Initializes WaterSnapshot class.

        Parameters
        ----------
        file : str
            Path to trajectory file.
        model : int
            Number of model in the trajectory.
        atom_lines : List[str]
            Water PDB lines.
        """
        super().__init__(file, model, atom_lines)
        self.water_coords, self.oxygen_coords = self.extract_water_coords(self.atom_lines)
        self.orientation = self.calculate_orientation(self.water_coords)
        self.voxel = None

    @staticmethod
    def extract_water_coords(water_lines):
        """
        Extracts coordinates from PDB lines.

        Parameters
        ----------
        water_lines : List[str]
            List of PDB water lines.

        Returns
        -------
            Water coordinates and oxygen coordinates.
        """
        water = [line[32:56].strip().split() for line in water_lines]
        return water, water[0]

    @staticmethod
    def calculate_orientation(water_coords):
        """
        Calculates orientation of the water molecules.

        Parameters
        ----------
        water_coords : List[str]
            Water coordinates extracted from PDB file.
        Returns
        -------
            Orientation vector of the water molecule.
        """
        coord_o, coord_h1, coord_h2 = water_coords
        coord_o = [float(c) for c in coord_o]
        coord_h1 = [float(c) for c in coord_h1]
        coord_h2 = [float(c) for c in coord_h2]

        v1 = np.subtract(coord_h1, coord_o)
        v2 = np.subtract(coord_h2, coord_o)

        v_middle = np.add(v1, v2)
        v_middle_magnitude = np.sqrt(
            v_middle[0] ** 2 + v_middle[1] ** 2 + v_middle[2] ** 2
        )

        orientation = [
            np.degrees(np.arccos(v_middle[0] / v_middle_magnitude)),
            np.degrees(np.arccos(v_middle[1] / v_middle_magnitude)),
            np.degrees(np.arccos(v_middle[2] / v_middle_magnitude)),
        ]

        return orientation

    def check_voxel(self, voxels):
        """
        Assigns the water snapshot to a specific voxel in the grid by measuring the distance of the oxygen coordinates
        to each voxel center.

        Parameters
        ----------
        voxels : List[Voxel]
            List of Voxel objects.
        """
        distances = {}

        for voxel in voxels:
            oxygen_coords = [float(c) for c in self.oxygen_coords]
            dist = np.subtract(oxygen_coords, voxel.voxel_center)
            dist = np.linalg.norm(dist)
            distances[voxel] = dist

        min_dist = min(list(distances.values()))
        result = [key for key, value in distances.items() if value == min_dist]

        self.voxel = result[0]


class Voxel:

    def __init__(self, voxel_center):
        """
        Initializes voxel class.
        Parameters
        ----------
        voxel_center : Tuple[float]
            Coordinates of voxel center.
        """
        self.voxel_center = voxel_center
        self.snapshots = []

    def add_snapshot(self, snap):
        """
        Adds snapshot to the Voxel.

        Parameters
        ----------
        snap : Snapshot
            Snapshot object containing information about file and model number.
        """
        self.snapshots.append(snap)

    def calculate_entropy(self):
        """
        Calculates entropy of water molecules in the Voxel.
        TODO: Alex
        """
        pass

    def calculate_enthalpy(self):
        """
        Calculates enthalpy of water molecules in the Voxel.
        TODO: Alex
        """
        pass


def extract_snapshots(simulation_output):
    """
    Extracts snapshots from the simulation output.

    Parameters
    ----------
    simulation_output : str
        Path to simulation output.

    Returns
    -------
        A list of WaterSnapshot objects.
    """
    print("Extracting water snapshots from {}.".format(simulation_output))
    all_files = glob.glob(os.path.join(simulation_output, "*", "trajectory*"))
    all_snapshots = []

    if len(all_files) == 0:
        raise FileNotFoundError(f"No output files detected in {simulation_output}.")

    for f in all_files:
        with open(f, "r") as file:
            file_content = file.read()
            matches = re.findall(r"(MODEL.+?ENDMDL)", file_content, re.DOTALL)

            for m in matches:
                model = re.search(r"MODEL\s+(\d+)", m).group(1)
                if model != "1":
                    water_match = re.findall(r".+HOH.+", m)
                    snapshot = WaterSnapshot(f, int(model), water_match)
                    all_snapshots.append(snapshot)

    return all_snapshots


def main(user_center, radius, output_simulation):
    grid = Grid(user_center, radius)
    grid.generate_voxels()
    snapshots = extract_snapshots(output_simulation)

    for snapshot in snapshots:
        snapshot.generate_properties()
        snapshot.check_voxel(grid.voxels)
        snapshot.voxel.add_snapshot(snapshot)  # Not sure if this makes sense
        grid.add_active_voxel(snapshot.voxel)
