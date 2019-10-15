import os
from string import Template
import argparse
import msm_pele.Helpers.helpers as hp
import msm_pele.Box.multiple_box as multiple_box
import msm_pele.Box.fixed_box as fixed_box
import msm_pele.Box.box_helpers as bhp
import msm_pele.constants as cs

__author__ = "Daniel Soler Viladrich"
__email__ = "daniel.soler@nostrumbiodiscovery.com"

# BOX CONSTANTS
KEYWORDS = ["MODEL", "RADIUS", "CENTER_X", "CENTER_Y", "CENTER_Z", "V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8"]
COORD = "{:>11.3f}{:>8.3f}{:>8.3f}"
CENTER = "{:.3f}"


def parseargs():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('bs', type=int, nargs='+', help='an integer for the accumulator')
    parser.add_argument('--points', '-p', type=str, help='File where to find all coord')
    args = parser.parse_args()
    return args.bs, args.points


def create_box(clusters, env, iteration):
    # Check exit simulation is properly finished
    sasa_column = "6" if env.noRMSD else "7"
    BS_sasa_min, BS_sasa_max = hp.is_exit_finish(env.adap_ex_output, env.test, criteria=sasa_column)
    # Retrieve BoxBuilderobj where
    # all abox carachteristics will be append
    box = BoxBuilder()
    if iteration == 0:
        # Retrieve info from pdb file
        if env.box:
            box.center, box.radius, box.box_type = box.from_file(env.box)
        # Retrieve info from user
        elif env.user_center and env.user_radius:
            box.center, box.radius, box.box_type = box.from_list(env.user_center, env.user_radius)
        # Create single box
        elif env.box_type == "fixed":
            box.center, box.radius, box.box_type = fixed_box.build_box(env.adap_ex_input, env.clusters_output, env.ligand_ref)
        # Create multiple box
        elif env.box_type == "multiple":
            box.center, box.radius, box.box_type = multiple_box.build_box(clusters, env)
        # Build box image
        box.to_pdb(env.box_temp)
        # Build box string for simulation
        box_string = box.to_PELE_string()
    else:
        env.box = os.path.join(env.pele_dir, "box.pdb")
        box.center, box.radius, box.box_type = box.from_file(env.box)
        # Build box string for simulation
        box_string = box.to_PELE_string()

    # Ensure box connectivity
    box.remove_clusters_out_of_box(env.clusters_output, file_path=env.cluster_output)
    return box_string, BS_sasa_min, BS_sasa_max


class BoxBuilder():

    def __init__(self):
        self.center = None
        self.radius = None
        self.box_type = None

    def from_file(self, input_file):
        """
        Retrieve info from MSM pdb box file
        """
        with open(input_file, 'r') as f:
            lines = hp.preproces_lines(f.readlines())

        try:
            center = [[float(line[5]), float(line[6]), float(line[7])] for line in lines if "CENTER" in line]
            radius = [float(line[2]) for line in lines if "RADIUS" in line]
        except ValueError:
            raise ValueError("{} not valid. Check the file is not a template.")

        box_type = "fixed" if len(radius) == 1 else "multiple"

        return center, radius, box_type

    def from_list(self, center, radius):
        """
        Retrieve Box Information
        from Python List
        """

        center = [center[i * 3: i * 3 + 3] for i in range(len(center))]
        radius = radius
        box_type = "multiple" if len(center) > 3 else "fixed"
        return center, radius, box_type

    def to_pdb(self, file):
        """
        Create a pdb file with the box charachteristics
        """
        # Initialize File
        with open(file, "w") as f:
            f.write("")

        for model, (center, radius) in enumerate(zip(self.center, self.radius), 1):
            # Build box Variables
            cx, cy, cz = center
            v1 = COORD.format(cx - radius, cy - radius, cz - radius)
            v2 = COORD.format(cx + radius, cy - radius, cz - radius)
            v3 = COORD.format(cx + radius, cy - radius, cz + radius)
            v4 = COORD.format(cx - radius, cy - radius, cz + radius)
            v5 = COORD.format(cx - radius, cy + radius, cz - radius)
            v6 = COORD.format(cx + radius, cy + radius, cz - radius)
            v7 = COORD.format(cx + radius, cy + radius, cz + radius)
            v8 = COORD.format(cx - radius, cy + radius, cz + radius)
            cx = CENTER.format(cx)
            cy = CENTER.format(cy)
            cz = CENTER.format(cz)
            values = [model, radius, cx, cy, cz, v1, v2, v3, v4, v5, v6, v7, v8]
            # Replace Template String
            replace = {keyword: value for keyword, value in zip(KEYWORDS, values)}
            box_tmp = Template(cs.BOX)
            box = box_tmp.safe_substitute(replace)
            # Write to a file
            with open(file, "a") as f:
                f.write(box)
            # Return File
        return file

    def to_PELE_string(self):
        """
        Build string for PELE control file
        """

        # Initialize Variables
        centers = self.center
        radiuses = self.radius
        box_type = self.box_type

        # Build String
        string = []
        if box_type == "fixed":
            string.append('"type": "sphericalBox",\n')
            for center, radius in zip(centers, radiuses):
                string.append('\n"radius": {},\n"fixedCenter":[{}]\n,'.format(radius, ",".join([str(coord) for coord in center])))
            string[-1] = string[-1].strip(",")
        elif box_type == "multiple":
            string.append('"type": "multiSphericalBox",\n"listOfSpheres":[')
            for center, radius in zip(centers, radiuses):
                string.append('{{\n"radius": {},\n"fixedCenter":[{}]\n}},'.format(radius, ",".join([str(coord) for coord in center])))
            string[-1] = string[-1].strip(",")
            string.append(']')
        return "\n".join(string)

    def remove_clusters_out_of_box(self, cluster_center, file_path=".", cluster_name="initial_{}.pdb"):
        """
        Remove clusters out of the exploration
        boxes. (They won't move on the simulation
        aggregating noise).
        """
        points = bhp.get_points(cluster_center)
        for i, point in enumerate(points):
            out_of_box = True
            for center, radius in zip(self.center, self.radius):
                cx, cy, cz = center
                x, y, z = point
                if ((x - cx)**2 + (y - cy)**2 + (z - cz)**2) <= (radius**2):
                    out_of_box = False
            if out_of_box:
                try:
                    os.remove(os.path.join(file_path, cluster_name.format(i)))
                except OSError:
                    pass
