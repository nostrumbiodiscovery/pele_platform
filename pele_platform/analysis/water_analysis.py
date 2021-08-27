from collections import defaultdict
import glob
import numpy as np
import os
import re
import math

from pele_platform.analysis.data import DataHandler
from pele_platform.Utilities.Helpers import helpers
from pele_platform.analysis.clustering import MeanShiftClustering
from pele_platform.analysis.analysis import Analysis


class WaterAnalysis:
    def __init__(
        self,
        simulation_output,
        be_column=4,
        trajectory="trajectory.pdb",
        report="report",
        topology=None,
        skip_initial_structures=True,
        water_ids=None,
    ):
        """
        Initializes WaterAnalysis class.

        Parameters
        -----------
        simulation_output : str
            Path to simulation output folder.
        be_column : int
            Column with energy metric, default 4.
        trajectory : str
            Trajectory name defaults to "trajectory.pdb", but you should use "trajectory.xtc" if using XTC format.
        report : str
            Report name.
        topology : str
            Path to the topology file (mandatory, if using XTC trajectories).
        skip_initial_structures : bool
            Toggle to skip extraction of the very first model.
        water_ids : list[tuple[str, int]]
            The list of water ids to track. Each water id is defined with
            a tuple that contains the PDB chain and the residue number
            corresponding to each water molecule to track.
        """
        self.output = simulation_output
        self.be_column = be_column
        self.trajectory = trajectory
        self.report = report
        self.topology = topology
        self.skip_initial_structures = skip_initial_structures
        self.water_ids = water_ids

        self.data_handler = DataHandler(
            sim_path=self.output,
            report_name=self.report,
            trajectory_name=self.trajectory,
            be_column=self.be_column,
            skip_initial_structures=self.skip_initial_structures,
        )

    @classmethod
    def from_parameters(cls, simulation_parameters):
        """
        Initializes WaterAnalysis class from simulation parameters.

        Parameters
        ----------
        simulation_parameters : Parameters
            Simulation parameters object.

        Returns
        --------
        WaterAnalysis object.
        """
        analysis = WaterAnalysis(
            simulation_output=os.path.join(
                simulation_parameters.pele_dir, simulation_parameters.output
            ),
            be_column=simulation_parameters.be_column,
            trajectory=simulation_parameters.traj_name,
            report=simulation_parameters.report_name,
            topology=simulation_parameters.topology,
            skip_initial_structures=not simulation_parameters.test,
            water_ids=simulation_parameters.water_ids_to_track,
        )

        return analysis

    def _write_waters(self, coordinates, clusters, path):
        """
        It writes the water molecules as a PDB file per each hydration site.
        Parameters
        ----------
        coordinates : a numpy.array with the coordinates of the atoms of the water molecules
        clusters: a numpy.array with a list of the cluster assigned to each water molecule
        path : string
            Output path of the directory where the molecules will be saved as a PDB file
        """
        nameatm = ["O  ", "H1 ", "H2 "]*round(len(coordinates)/3)
        count = 0
        nclusters = np.amax(clusters) + 1

        if not os.path.exists(path):
            os.makedirs(path)

        for i in range(nclusters):
            waters_id = np.where(clusters == i)[0]
            name_f = "water_molecules_{:03d}.pdb".format(i)
            path_f = os.path.join(path,name_f)
            
            with open(path_f, 'w+') as f:
                for label in waters_id:
                    label = int(label)
                    for j in range(3):
                        atompos = coordinates[3*label + j]
                        f.write("ATOM    {:3d}  ".format(3*label + j) + nameatm[3*label + j] +
                            " HOH A {:3d} ".format(label) +
                            "{:>11.3f}{:>8.3f}{:>8.3f}  ".format(*atompos) +
                            "1.00 1.00\n")

    def quaternion_rotation_vector_3d(self,v, q):
        """
        This function calculates the result of the rotation of a 3d vector by a quaternion

        Parameters
        ----------
        v: np.array of dimension 3
            3d position vector with componentns x,y,z
        q: np.array of dimension 4
            quaternion representing an orientation

        Returns
        -------
        v_rot: np.array of dimension 3
            3d position vector resulting from the rotation of v by q
        """
        x0 = v[0]; y0 = v[1]; z0 = v[2]; w = q[0]
        x = q[1]; y = q[2]; z = q[3]

        q_aux = np.array([x0*x + y0*y + z0*z,
                          x0*w - y0*z + z0*y,
                          x0*z + y0*w - z0*x,
                         -x0*y + y0*x + z0*w])

        w2 = q_aux[0]; x2 = q_aux[1]; y2 = q_aux[2]; z2 = q_aux[3]

        v_rot = np.array([w*x2 + x*w2 + y*z2 - z*y2,
                          w*y2 - x*z2 + y*w2 + z*x2,
                          w*z2 + x*y2 - y*x2 + z*w2])
        return v_rot

    def quaternion_product(self,q1,q2):
        """
        This function calculates the product of two quaternions, combining the two rotations

        Parameters
        ----------
        q1, q2: np.arrays of dimension 4
                 quaternions representing two rotations

        Returns
        -------
        qT: np.array of dimension 4
             quaternion representing the total rotation
        """
        w1 = q1[0]; x1 = q1[1]; y1 = q1[2]; z1 = q1[3]
        w2 = q2[0]; x2 = q2[1]; y2 = q2[2]; z2 = q2[3]

        qT = np.array([w2*w1 - x2*x1 - y2*y1 - z2*z1,
                      w2*x1 + x2*w1 + y2*z1 - z2*y1,
                      w2*y1 - x2*z1 + y2*w1 + z2*x1,
                      w2*z1 + x2*y1 - y2*x1 + z2*w1])

        return qT

    def conjugate_quaternion(self,q):
        """
        This function calculates the conjugate of a quaternion

        Parameters
        ----------
        q: np.array of dimension 4

        Returns
        -------
        qc: np.array of dimension 4
        """
        w = q[0]; x = q[1]; y = q[2]; z = q[3]
        qc = np.array([w,-x,-y,-z])
        return qc

    def water_mass_centre(self,position):
        """"
        This function calculates the mass centre of a water molecule

        Parameters
        ----------
        position: np.array of dimension (3,3)
                   positions (x,y,z) of atoms O, H, H

        Returns
        -------
        mass_centre: np.array of dimension 3
        """

        relative_mass_O = 15.999/18.015; relative_mass_H = 1.008/18.015

        pos_O = position[0]
        pos_H1 = position[1]
        pos_H2 = position[2]

        mass_centre = relative_mass_O*pos_O + relative_mass_H*(pos_H1 + pos_H2)

        return mass_centre

    def water_orientation_quaternion(self,position):
        """
        This function calculates the quaternion associated to the orientation of a water molecule

        Parameters
        ----------
        position: np.array of dimension (3,3)
                   coordinates (x,y,z) of atoms O, H, H

        Returns
        -------
        qT: np.array of dimension 4
             quaternion associated to the orientation of the molecule
        """
        relative_mass_O = 15.999/18.015; relative_mass_H = 1.008/18.015

        pos_O = position[0]
        pos_H1 = position[1]
        pos_H2 = position[2]

        mass_centre = relative_mass_O*pos_O + relative_mass_H*(pos_H1 + pos_H2)

        vec1 = pos_O - mass_centre
        vec1 = vec1 / np.linalg.norm(vec1)

        vec2 = pos_H1 - mass_centre
        vec2 = vec2 / np.linalg.norm(vec2)

        # Parameters to calculate the orientation of the water molecule
        phi1 = np.arccos(vec1[0])
        d1 = math.sqrt(vec1[1]**2 + vec1[2]**2)
        q1 = np.array([np.cos(phi1/2), 0, (vec1[2]/d1)*np.sin(phi1/2), -(vec1[1]/d1)*np.sin(phi1/2)]) #First quaternion

        vec2_r = self.quaternion_rotation_vector_3d(vec2,q1)
        vec2_r = vec2_r / np.linalg.norm(vec2_r)

        d2 = math.sqrt(vec2_r[1]**2 + vec2_r[2]**2)
        phi2 = np.arccos(vec2_r[1]/d2)

        q2 = np.array([np.cos(phi2/2),np.sin(phi2/2),0,0])
        qT = self.quaternion_product(q2,q1)
        qT = self.conjugate_quaternion(qT)

        return qT

    def entropy_clusters(self,coordinates, clusters, watersites_data):

        """
        This function calculates the orientational entropy associated to each hydration site

        Parameters
        ----------
        coordinates: numpy.array of [M, N, D] shape
                        M - number of models
                        N - number of atoms in the molecule (3)
                            Obs: It will actually be 3*nwaters
                        D - number of dimensions (3)
        clusters: numpy.array
                   labels of the cluster assigned to each water molecule
        watersites_data: List[dictionary]
                          information of each hydration site, contained in a dictionary
                          "cluster"(int), "position"(np.array (3,3) ), "population" (float)

        Returns
        -------
        clusters_entropy: numpy.array of len(watersites_data) with the entropy of each hydration site
        """

        #Constants
        models, atoms, dimensions = coordinates.shape
        nwaters = round(atoms/3)
        rho = 0.0327594 #Number density molecules/A^3 (A of Angstrom) of bulk water - depends on the used model
        Rgas = 1.9858775*1e-3 #Gas constant in kcal/mol·K, it converts the entropy to thermodynamic units
        relative_mass_O = 15.999/18.015; relative_mass_H = 1.008/18.015

        coordinates = coordinates.reshape(models*atoms, dimensions)
        coordinates = coordinates.reshape(round(models*atoms/3),3, dimensions)

        water_distances = list()

        ref_orientations = [None]*len(watersites_data)

        nclusters = len(watersites_data)
        clusters = clusters.astype(int)

        # Computing the reference water orientation
        for i in range(nclusters):
            molecules = np.where(clusters==i)[0] #Look for the water molecules in this cluster
            hydration_site = next(item for item in watersites_data if item["cluster"] == i) #Find data from this cluster
            hydration_centre = hydration_site["position"]
            position_ref = np.zeros((3,3)) # Find the average position of H, O, O
            for x in molecules:
                position_ref = position_ref + coordinates[x]
            position_ref = position_ref / len(molecules)
            ref_orientations[i] = self.water_orientation_quaternion(position_ref)


        water_distances = list() #List where we will save the entropy metric "distance" of each water molecule

        for i,position in enumerate(coordinates):
            cluster = clusters[i]
            q = self.water_orientation_quaternion(position)
            q = q / np.linalg.norm(q) # Normalize quaternion again for precision
            
            hydration_site = next(item for item in watersites_data if item["cluster"] == cluster)
            hydration_centre = hydration_site["position"]

            q_ref = ref_orientations[cluster]
            q_ref = q_ref / np.linalg.norm(q_ref)

            mass_centre = self.water_mass_centre(position)

            # Translational distance between molecule and hydration centre
            d_trans = np.linalg.norm(hydration_centre - mass_centre)

            dotprod = np.dot(q,q_ref)

            # Fix for numerical precision
            if(dotprod > 1):
                dotprod = 1.0
            elif(dotprod < -1):
                dotprod = -1.0

            # Orientational distance between molecule and reference
            d_orient = 2*np.arccos(np.absolute(dotprod))

            # Total distance
            d_tot = math.sqrt(d_trans**2 + d_orient**2)

            water_distances.append(d_tot)

        water_distances = np.asarray(water_distances)
        #At this point we should have an array water_distances with the same length as the variable 'clusters'

        clusters_entropy = np.zeros(nclusters)

        # Calculate the entropy for each cluster
        for i in range(nclusters):
            molecules = np.where(clusters==i)[0] # Molecules of this cluster
            n = len(molecules)
            sum_contributions = 0
            for x in molecules:
                sum_contributions = sum_contributions + np.log((water_distances[x]**6)*n*np.pi*rho/48)    
            npairs = 1
            pair_contributions = 0
            if (n>1):
                npairs = n*(n-1)/2
                for j in range(n):
                    for k in range(j+1,n):
                        dpair = water_distances[molecules[j]]**2 + water_distances[molecules[k]]**2
                        pair_contributions = pair_contributions + np.log(npairs*(dpair**6)*((np.pi*rho)**2)/46080)
            clusters_entropy[i] = Rgas*(sum_contributions/n) + Rgas*(pair_contributions/npairs)

        # Change to energetic units
        clusters_entropy = -300*clusters_entropy

        return clusters_entropy

    def _write_entropy(self, entropies, path):
            """
            It writes the entropy of each hydration site in a txt file

            Parameters
            ----------
            entropies: a numpy.array with the entropy associated to each hydration site
            path : string
                Output path of the directory where the molecules will be saved as a PDB file
            """

            if not os.path.exists(path):
                os.makedirs(path)

            name_f = "clusters_entropy.txt"
            path_f = os.path.join(path,name_f)

            with open(path_f, 'w+') as f:
                f.write("Entropic contributions -T·dS kcal/mol] \n")
                for count, value in enumerate(entropies):
                    f.write("CLUSTER    {:3d}:  ".format(count) + "{:>11.3f}".format(value) + "\n")  

 
    def analyse_waters(self, path):
        """
        Runs the whole water analysis workflow:
            - Creates folder.
            - Extracts water coordinates (only tracked IDs).
            - Clusters using mean shift algorithm.
            - TODO: Analyses energy and entropy of each water site.
        """
        helpers.check_make_folder(path)

        coordinates = self.extract_water_coordinates(
            skip_initial_structures=self.skip_initial_structures
        )

        models, atoms, dimensions = coordinates.shape


        clusters, water_sites = self.get_water_sites(coordinates, path)
        print(f"Generated {len(clusters)} clusters.")

        waterclusters = np.concatenate((np.array([np.arange(len(clusters))]).T, np.array([clusters]).T), axis=1)
        waterclusters = np.repeat(waterclusters, 3, axis=0)

        with open(os.path.join(path, "water_coords.txt"), "w+") as f:
            np.savetxt(f, np.concatenate((waterclusters,coordinates.reshape(models*atoms, dimensions)), axis=1))
        
        self._write_waters(coordinates.reshape(models*atoms, dimensions),clusters,os.path.join(path,"water_molecules/"))
        
        clusters_entropy = self.entropy_clusters(coordinates, clusters, water_sites)
        self._write_entropy(clusters_entropy, path)

    def extract_water_coordinates(self, skip_initial_structures=True):
        """
        Extracts water coordinates from all trajectories.

        Returns
        -------
        """
        trajectories = glob.glob(
            os.path.join(self.output, "*", self.trajectory.replace(".", "*"))
        )

        if self.topology:
            coordinates = self.extract_waters_from_XTC(
                trajectories,
                water_ids=self.water_ids,
                skip_initial_structures=skip_initial_structures,
            )
        else:
            coordinates = self.extract_waters_from_PDB(
                trajectories,
                water_ids=self.water_ids,
                skip_initial_structures=skip_initial_structures,
            )
        print("Water coordinates extracted.")
        return coordinates

    @staticmethod
    def extract_waters_from_PDB(all_files, water_ids, skip_initial_structures=True):
        """
        Extracts coordinates of water molecules from PDB trajectories.

        Parameters
        ----------
        all_files : List[str]
            List of all trajectory files.
        skip_initial_structures : bool
            Toggle to exclude the first model in each trajectory.
        water_ids : list[tuple[str, int]]
            The list of water ids to track. Each water id is defined with
            a tuple that contains the PDB chain and the residue number
            corresponding to each water molecule to track.

        Returns
        -------
            An array of [M, N, D] shape, where M - number of models, N - number of atoms in the molecule (3), D - number
            of dimensions (3).
        """
        water_coordinates = defaultdict(list)

        if len(all_files) == 0:
            raise FileNotFoundError("No output files detected.")

        for f in all_files:
            with open(f, "r") as file:
                file_content = file.read()
                models = re.findall(r"(MODEL.+?ENDMDL)", file_content, re.DOTALL)

                for model in models:
                    model_number = re.search(r"MODEL\s+(\d+)", model).group(1)
                    if model_number == 1 and skip_initial_structures:
                        continue

                    water_lines = re.findall(r".+HOH.+", model)

                    for line in water_lines:
                        chain = line[21]
                        residue_number = line[22:26].strip()

                        if (chain, int(residue_number)) in water_ids:
                            try:
                                x, y, z = (
                                    float(line[30:38]),
                                    float(line[38:46]),
                                    float(line[46:54]),
                                )
                            except ValueError:
                                print(f"Warning: invalid PDB format in trajectory {f}.")

                        point = np.array((x, y, z))
                        water_coordinates[
                            f"{f}_{model_number}_{chain}_{residue_number}"
                        ].append(point)

        all_coordinates = [
            np.array(coordinates) for coordinates in water_coordinates.values()
        ]
        all_coordinates = np.array(all_coordinates)

        n_models, n_atoms, dimensions = all_coordinates.shape
        assert dimensions == 3, "Extracted water coordinates are not 3-dimensional."
        assert (
            n_atoms == 3
        ), "Water molecules do not seem to have the correct number of atoms."

        return all_coordinates

    @staticmethod
    def extract_waters_from_XTC(all_files, water_ids, skip_initial_structures=True):
        """
        Extracts coordinates of water molecules from XTC trajectories.

        Parameters
        ----------
        all_files : List[str]
            List of all trajectory files.
        skip_initial_structures : bool
            Toggle to exclude the first model in each trajectory.
        water_ids : list[tuple[str, int]]
            The list of water ids to track. Each water id is defined with
            a tuple that contains the PDB chain and the residue number
            corresponding to each water molecule to track.

        Returns
        -------
            An array of [M, N, D] shape, where M - number of models, N - number of atoms in the molecule (3), D - number
            of dimensions (3).
        """
        water_coordinates = list()

        # TODO: Implement extraction from XTC trajectories, preferably as a static method.

        return water_coordinates

    def filter_oxygen_coordinates(self, coordinates):
        n_models, n_atoms, dimensions = coordinates.shape
        oxygens = round(n_atoms/3)
        new_coordinates = np.zeros((n_models,oxygens,dimensions))
        for i in range(n_models):
            for j in range(oxygens):
                new_coordinates[i][j] = coordinates[i][3*j]

        return new_coordinates

    @staticmethod
    def get_water_sites(coordinates, path):
        """
        Cluster water coordinates using mean shift algorithm and write the centroids to PDB file.

        Parameters
        ----------
        coordinates : List[array]
            List of numpy arrays with water coordinates.
        path : str
            Path to save the PDB files.
        """
        watersites_data = list()

        clustering = MeanShiftClustering(2.0)  # hard=coded bandwidth value for water
        
        n_models, n_atoms, dimensions = coordinates.shape
        oxygens = round(n_atoms/3)
        new_coordinates = np.zeros((n_models,oxygens,dimensions))
        for i in range(n_models):
            for j in range(oxygens):
                new_coordinates[i][j] = coordinates[i][3*j]

        water_clusters, estimator = clustering.get_clusters(new_coordinates)

        populations = Analysis._get_cluster_populations(water_clusters)
        output_path = os.path.join(path, "water_sites.pdb")
        Analysis._write_centroids(populations, estimator, output_path)
        centroids = estimator.cluster_centers_

        for cluster, centroid in enumerate(centroids):
            watersites_data.append({
                "cluster": cluster,
                "position": centroid,
                "population": populations[cluster]
            })

        return water_clusters, watersites_data
