Input flags documentation
###########################

Compulsory flags PELE
--------------------------------

- **system**: Path to the input pdb file cointaining ligand and receptor in your desired initial conformation (except for a global exploration)

 
- **residue**: Residue name of the ligand to be perturbed


- **chain**: Chain of the ligand to be perturbed


- **cpus**: Cpus to use

..  code-block:: yaml

    system: "/home/daniel/PR_complex.pdb"
    residue: "LIG"
    chain: "L"
    cpus: 200


Compulsory flags FragPele
-----------------------------

Frag PELE grows an atom onto a core in N growing steps while moving protein and ligand.
Afterwards a final sampling simulation is run to fully explore the ligand-protein conformational space.

- **frag_core**: Core of the molecule we want to add fragments to. Required parameter

- **Method to use**: Choose on of the available methos. For more please refer here.

- **resname**: Residue name of the frag_core ligand

- **cpus**: Cpus to use. Default=48

..  code-block:: yaml

    frag_core: "/home/daniel/PR_core.pdb"
    frag_input: "/home/daniel/serie_file.conf"
    frag_core: "LIG"
    cpus: 48


Optative flags
----------------------------


General settings
====================

Configure the settings of the simulation and the path to all dependencies in case of need (non-default installation).

- **test**: Run a quick test to check the simulation works (~2 min). **Never use the control files from the test as input for a production simulation as temperature, ANM and minimization are twicked to made the simulation faster!!!!**
 
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

- **log**: Retrieve PELE logfiles during simulation. Default=False

- **verbose**: Set to true to activate verbose mode in PELE. DEfault=False

- **anm_freq**: Every how many steps to perform anm. Default=4

- **anm_displacement**: Angstrom to displace carbon alphas in each ANM movement. Default=0.75

- **anm_modes_change**: Number of steps before we change to a new normal mode movement. Default=4

- **sidechain_freq**: Every how many steps to perform sidechain sampling. Default=2

- **min_freq**: Every how many steps to perform minimization. Default=1

- **water_freq**: Every how many steps to perform water perturbation. Default=1

- **temperature**: Temperature of the simulation. Default=1500

- **solvent**: Solvent of the simulation. (OBC or VDGBNP). Default=VDGBNP

- **sidechain_res**: Receptor sidechain resolution. Default=10

- **overlap_factor**: Vanderwals overlap factor (More in PELE docs). Default=0.65

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

- **debug**: Use this flag to only create the inputs of the simulation. No simulation is run. (Usefull to transport it to another machine)

- **spawning**: Spawning type ([independent, inverselyProportional or epsilon so far]). Default: inverselyProportional

- **density**: Density type ([null, exitContinuous...]. More in AdaptivePELE docs). Default: null

- **cluster_values**: Clusterization values. More in AdaptivePELE. Default: Depending on simulation type

- **cluster_conditions**: Clusterization condition. More in AdaptivePELE. Default: Depending on simulation type

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

FragPELE params
===================

**These flags are exclusive of the FragPele modes not PELE**

- **growing_steps**: Number of steps to grow the fragment with.

- **steps_in_gs**: Number of pele steps within each growing step

- **sampling_steps**: Number of pele steps in the final sampling simulation

- **protocol**: Type of protocol. options = [HT, ES]. For more info please refere here.


..  code-block:: yaml

    growing_steps: 6
    steps_in_gs: 6
    sampling_steps: 20
    protocol: HT
    cpus: 24

PPI params
===============

**These flags are exclusive of the ppi: true mode**

- n_components: Number of clusters after global exploration. In other words, number of inputs for the refinment exploration after the global simulation. Default: 10


..  code-block:: yaml

    n_components: 10


Constraints
==================

This section allows the user to change the constraint values.

- **ca_constr**: Carbon alpha constraints. Default=0.5

- **interval_constr**: Every how many carbon alphas to apply the constraints. Default:10

- **metal_constr**: Metal constraints. Default=200

- **water_constr**: Water constraints. Default=5

- **constrain_smiles**: SMILES string to indicate what part of the molecule to constraint. Default=None

- **smiles_constr**: Numeric value of the SMILES constraints. Default=10

- **external_constraints**: You can specify 2 types of constraints. Positional constraints or atom-atom constraint. (Example below)

  - The positional constraints are given either by: 
        - springConstant-atomnumber. i.e. "10-17"
        - springConstant-chain:resnum:atomname. i.e. "5-A:1:H"

  - The atom-atom constraints are specified either by: 
        - springConstant-equilibriumDistance-atomnumber1-atomnumber2. i.e. "50-2.34-17-4159"
        - springConstant-equilibriumDistance-chain1:resnum1:atomname1-chain2:resnum2:atomname2. i.e. "50-2.34-A:1:H-L:1:C21"

..  code-block:: yaml

    ca_constr: 2
    interval_constr: 10
    metal_constr: 100
    water_constr: 5
    constrain_smiles: "C2CCC1CCCCC1C2"
    smiles_constr: 5
    external_constraints:
    - "10-17" #constraint of 10kcal/mol at atomnumber 17
    - "5-A:1:H" ##constraint of 10kcal/mol at atom with chain A residuenumber 1 and atomname H
    - "50-2.34-17-4159" #constraint of 50kcal/mol with equilibrium distance of 2.34 between atomnumbers 17 & 4159
    - "50-2.34-A:1:H-L:1:C21" #constraint of 50kcal/mol with equilibrium distance of 2.34 between atoms with respective chain resnum and atomname


WaterPerturbation
======================


    - **water_exp**: Exploration of the hydratation sites of a binding site by perturbing and clusterizing a single water. More advance features will be later implemented to discriminate between "happy" and "unhappy" waters.

    - **water_lig**: Perturb one or several water molecules while exploring the conformational space of the ligand.

Example water exploration:

..  code-block:: yaml

  water_exp:
    - M:1
    - M:2

Example water ligand:

..  code-block:: yaml

    water_lig:
    - M:1
    - M:2

Simulation Parameters
========================

- **box_water**: Center of the box for the waters. Default: Centroid of the center of masses of all water molecules.

- **water_radius**: Radius of the water box. Default=7

- **water_trials**: Numerical trials on water perturbation. Default=10000

- **water_constr**: COM constraint applied to th water molecule after perturbation. Default=0

- **water_temp**: Temperature of the water perturbation step. Default=5000

- **water_overlap**: Overlap factor of water. Default=0.78


..  code-block:: yaml

    box_water:
    - 20
    - 30
    - 20
    water_radius: 8
    water_trials: 500
    water_constr: 0.5
    water_temp: 2000
    water_overlap: 0.5


Metrics
=============

Metrics to track along the simulation

- **atom_dist**: Calculate distance between two atomnumbers. To calculate more than one append them in column as the example below. Default=None

    - The atomdist can be specified via chain:resnum:atomname i.e. A:2:CA

- **rmsd_pdb**: Calculate rmsd of the ligand to a native pdb structure


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

- **mae**: To extract the best energy and cluster poses as .mae files with the metrics as properties (schrodinger need it). Default: false

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
