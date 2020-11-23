All PELE packages
#####################################

General settings
====================

Configure the settings of the simulation and the path to all dependencies in case of need (non-default installation).

- **test**: Run a quick test to check the simulation works (~2 min). **Never use the control files from the test as input for a production simulation as temperature, ANM and minimization are twicked to make simulation faster!!!!**

- **usesrun**: Use srun binary to run PELE. Only when using intel processors.

- **pele_exec**: Use a pele executable that is not the default one. **Needs to be used with pele_data and pele_documents**. default: $PELE/bin/Pele_mpi

- **pele_data**: Use a pele data folder that is not the default one. default: $PELE/Data

- **pele_documents**: Use a pele documents folder that is not the default one. default: $PELE/Documents

- **pele_license**: Use a pele_license path that is not the default one. default: $PELE/licenses

- **schrodinger**: Use a schrodinger path that is not the default one. default: $SCHRODINGER

..  code-block:: yaml


  test: true
  usesrun: false
  pele_exec: "/home/pele/bin/Pele_mpi"
  pele_data: "/home/pele/Data/"
  pele_documents: "/home/pele/Documents/"
  pele_license: "/home/pele/licenses"
  schrodinger: "/home/pele/schrodinger2020-1/"


Receptor preparation
=======================

Configure the parameters of the PPP (Protein Pele Preparation)

- **skip_preprocess**: Skip protein pele preparation. Default: False

- **noTERs**: Don't include TERs on preparation. Used if PPP gets confuse with insertion codes or other. Default: False

- **charge_ters**: Charge terminals of the protein. Default: False

- **nonstandard**: List of names of nonstandard residues that will be omitted in protein pele preparation. Default=[]

- **prepwizard**: Run Prepwizard (Still on testing version). Default: False

..  code-block:: yaml

  preprocess_receptor: true
  noTERs: false
  charge_ters: false
  nonstandard:
    - TPO
  prepwizard: false


Ligand preparation
======================

Configure the parameters of the PlopRotTemp to extract the ligand forcefield parameters.

- **gridres**: Resolution of the rotamers when sampling. Default: 10 degrees

- **core**: Atomnumber of the atom that will be included as part of the rigid core. Default=None

- **maxtorsion**: Maximum number of rotamers per flexible sidechain. Default: 4

- **n**: Maximum number of flexible sidechains in a molecule, Default: None

- **mae_lig**: Mae file to extract the cuantum charges from. Default: None

- **template**: External forcefield templaters

- **rotamers**: External rotamer libraries

- **skip_ligand_prep**: Skip preparation of that resiude. This could be usefull to bypass problems with PlopRotTemp when creating the ligand parameters.


..  code-block:: yaml

  gridres: 10
  core: -1
  maxtorsion: 4
  n: 5
  mae_lig: "/home/dsoler/lig.mae"
  templates:
    - "/home/dsoler/mgz"
    - "/home/dsoler/ligz"
  rotamers:
    - "/home/dsoler/MG.rot.assign"
    - "/home/dsoler/LIG.rot.assign"
  skip_ligand_prep:
    - "LIG"

Box parameters
=================

Parameters to set the exploration Box:

- **box_radius**: Radius of the box. Default=[induced_fit (10), local_exploration (30), global_exploration (50)]

- **box_center**: Center of the box. Default=[indeuced_fit&local_exploration (CM of the ligand), global (calculater center)]


..  code-block:: yaml

  box_radius: 30
  box_center:
    - 20
    - 30
    - 50


Simulation params
====================

- **seed**: Seed of the job for reproducibility. Default=12345

- **log**: Retrieve PELE log files during simulation. Default=False

- **verbose**: Set to true to activate verbose mode in PELE. Default=False

- **anm_freq**: Every how many steps to perform anm. Default=4

- **anm_displacement**: Angstrom to displace carbon alphas in each ANM movement. Default=0.75

- **anm_modes_change**: Number of steps before we change to a new normal mode movement. Default=4

- **sidechain_freq**: Every how many steps to perform sidechain sampling. Default=2

- **min_freq**: Every how many steps to perform minimization. Default=1

- **water_freq**: Every how many steps to perform water perturbation. Default=1

- **temperature**: Temperature of the simulation. Default=1500

- **solvent**: Solvent of the simulation. (OBC or VDGBNP). Default=VDGBNP

- **sidechain_res**: Receptor sidechain resolution. Default=10

- **overlap_factor**: van der Wals overlap factor (More in PELE docs). Default=0.65

- **steric_trials**: Number of steric trials (More in PELE docs). Default=250

..  code-block:: yaml

  seed: 312312
  log: true
  verbose: true
  anm_freq: 4
  anm_displacement: 0.5
  anm_modes_change: 3
  sidechain_freq: 2
  min_freq: 1
  water_freq: 1
  temperature: 1500
  solvent: "VDGBNP"
  sidechain_res: 30
  overlap_factor: 0.65
  steric_trials: 250



PELE params
===================

**These flags are exclusive of the PELE modes not fragPELE**

- **iterations**: Adaptive epochs to run. Set to 1 by default if using PELE

- **steps**: Pele steps in each iteration

- **debug**: Use this flag to only create the inputs of the simulation. No simulation is run. (Useful to transport it to another machine)

- **spawning**: Spawning type ([independent, inverselyProportional or epsilon so far]). Default: inverselyProportional

- **density**: Density type ([null, exitContinuous...]. More in AdaptivePELE docs). Default: null

- **cluster_values**: Clustering values. More in AdaptivePELE. Default: Depending on simulation type

- **cluster_conditions**: Clustering condition. More in AdaptivePELE. Default: Depending on simulation type

- **equilibration**: Whether to run initial equilibration or not. Default: false

- **equilibration_steps**: Equilibration steps. Default: 2

- **adaptive_restart**: Use adaptive restart with the working folder option to restart the simulation. Default: false

- **report**: Change the name of the report file. Default: report

- **traj**: Change the name of the trajectory file. Default: trajectory.pdb

..  code-block:: yaml

    iterations: 30
    steps: 12
    debug: true
    spawning: "epsilon"
    density: "exitContinuous"
    cluster_values: [2,3,4]
    cluster_conditions: [0.8, 0.6, 0.2]
    equilibration: false
    equilibration_steps: 10
    adaptive_restart: true
    working_folder: "folder_to_restart"
    report: report
    traj: trajectory.xtc


Constraints
==================

This section allows the user to change the constraint values.

- **ca_constr**: Carbon alpha constraints. Default=0.5

- **interval_constr**: Every how many carbon alphas to apply the constraints. Default:10

- **water_constr**: Water constraints. Default=5

- **constrain_smiles**: SMILES string to indicate what part of the molecule to constrain. Default=None

- **smiles_constr**: Numeric value of the SMILES constraints. Default=10

- **external_constraints**: You can specify 2 types of constraints. Positional constraints or atom-atom constraint. (Example below)

  - The positional constraints are given either by:
        - springConstant-atomnumber. i.e. "10-17"
        - springConstant-chain:resnum:atomname. i.e. "5-A:1:H"

  - The atom-atom constraints are specified either by:
        - springConstant-equilibriumDistance-atomnumber1-atomnumber2. i.e. "50-2.34-17-4159"
        - springConstant-equilibriumDistance-chain1:resnum1:atomname1-chain2:resnum2:atomname2. i.e. "50-2.34-A:1:H-L:1:C21"

- **remove_constraints**: Do not place constraints on the carbon-alpha of the protein. Default: False


..  code-block:: yaml

    ca_constr: 2
    interval_constr: 10
    water_constr: 5
    constrain_smiles: "C2CCC1CCCCC1C2"
    smiles_constr: 5
    external_constraints:
    - "10-17" #constrain of 10 kcal/mol at atom number 17
    - "5-A:1:H" ##constrain of 10 kcal/mol at atom with chain A residue number 1 and atomname H
    - "50-2.34-17-4159" #constrain of 50 kcal/mol with equilibrium distance of 2.34 between atom numbers 17 & 4159
    - "50-2.34-A:1:H-L:1:C21" #constrain of 50 kcal/mol with equilibrium distance of 2.34 between atoms with respective chain resnum and atomname
    remove_constraints: true

Metal constraints
+++++++++++++++++++++

Algorithm to automatically set metal constraints around the ligand.

- **no_metal_constraints**: Ignore all metals in the PDB file, no constraints will be set automatically. Default=False

- **permissive_metal_constr**: Expand the search for coordinated atoms by allowing 35% deviation from “ideal” angles. If the algorithm finds a valid geometry it will include the metal constraint into the simulation. Default=False

- **constrain_all_metals**: Constrain all atoms around the metal, regardless of the angles or coordination number. Default=False

- **external_constraints**: Set a manual constraint containing a metal atom to disable search for this particular metal. Default=[]


..  code-block:: yaml

    no_metal_constraints: true
    permissive_metal_constr: true
    constrain_all_metals: true
    external_constraints:
        - "50-2.34-A:1:H-L:1:MG" #constrain of 50kcal/mol with equilibrium distance of 2.34 between atoms with respective chain resnum and atomname

Metal polarisation
====================

An optional flag to adjust charges on the metals by dividing them by certain factor.

- **polarize_metals** - adjust charges on the metals by dividing them by 2 (unless other value is set in polarization_factor)

- **polarization_factor** - factor by which the metal charges should be divided

..  code-block:: yaml

    polarize_metals: true
    polarization_factor: 2 # Mg2+ will have a charge of +1

Metrics
=============

Metrics to track along the simulation

- **atom_dist**: Calculate distance between two atom numbers. To calculate more than one append them in column as the example below. Default=None

    - The atom distances can be specified via chain:resnum:atomname, e.g. A:2:CA.

- **rmsd_pdb**: Calculate RMSD of the ligand to a native pdb structure


..  code-block:: yaml

    atom_dist:
        # Distance between the A:2:CA and B:3:CG also between A:5:N and B:3:CG. Append more if desired.
        - "A:2:CA"
        - "B:3:CG"
        - "A:5:N"
        - "B:3:CG"
    rmsd_pdb: "/home/dsoler/native.pdb"


Analysis
=============

Run a post simulation analysis to extract plots, top poses and clusters.

- **only_analysis**: Analyse PELE simulation without running it.

- **analysis_nclust**: Numbers of clusters out of the simulation. Default: 10

- **be_column**: Column of the binding energy in the reports starting by 1. Default: 5

- **te_column**: Column of the total energy in the reports starting by 1. Default: 4

- **limit_column**: Specify the column where your external metrics start. Default: 6

- **mae**: To extract the best energy and cluster poses as .mae files with the metrics as properties (Schrodinger needs it). Default: false

- **analysis**: Whether to run or not the analysis at the end of the simulation. Default: true

..  code-block:: yaml

    only_analysis: true
    be_column: 5
    te_column: 4
    limit_column: 6
    mae: true

Output
==========

Configure the output

- **working_folder**: Name of the main working folder where to store the processed input, control files and the simulation folder. Default="resname_Pele_X" where X is a number.

- **output**: Output folder of the simulation. Default=output

..  code-block:: yaml

    working_folder: "NOR_solvent_OBC"
    output: "output_sim"
