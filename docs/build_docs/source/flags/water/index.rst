Water perturbation
===================

- **n_waters**: Number of waters to randomly add into your simulation and perturb. Default=0

- **waters**: Water molecules to be perturbed in AquaPELE steps. Users can indicate specific water IDs, e.g. "W:15" or select "all_waters" option to perturb all water molecules present in the system.

- **box_water**: Center of the box for the waters. Default: Centroid of the center of masses of all water molecules.

- **water_radius**: Radius of the water box. Default=7

- **water_trials**: Numerical trials on water perturbation. Default=10000

- **water_constr**: COM constrain applied to th water molecule after perturbation. Default=0

- **water_temp**: Temperature of the water perturbation step. Default=5000

- **water_overlap**: Overlap factor of water. Default=0.78


..  code-block:: yaml

    n_waters: 3 # Compulsory, if no water molecules are present in the system
    box_water:
    - 20
    - 30
    - 20
    water_radius: 8
    water_trials: 500
    water_constr: 0.5
    water_temp: 2000
    water_overlap: 0.5
    # waters: "all_waters" # to perturb all waters in the system
    # waters:
        - "W:15" # chain ID and residue number
        - "W:21"
