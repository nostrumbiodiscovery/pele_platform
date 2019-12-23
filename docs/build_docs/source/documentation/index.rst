Input Parameters
######################

Compulsory flags
--------------------

- **system**: Path to the input pdb file cointaining ligand and receptor in your desired initial conformation (except for a global exploration)

 
- **residue**: Residue name of the ligand to be perturbed


- **chain**: Chain of the ligand to be perturbed


- **cpus**: Cpus to use

..  code-block:: yaml

    system: "/home/daniel/PR_complex.pdb"
    residue: "LIG"
    chain: "L"
    cpus: 200


Optative flags
-------------------

Job parameters
=================

Configure the main important parameters for the job


- **iterations**: Adaptive epochs to run. Set to 1 by default if using PELE

- **steps**: Pele steps in each iteration

- **test**: Run a quick test to check the simulation works (~2 min). **Never use the control files from the test as input for a production simulation as temperature, ANM and minimization are twicked to made the simulation faster!!!!**
 
- **usesrun**: Use srun binary to run PELE. Only when using intel processors.

- **debug**: Use this flag to only create the inputs of the simulation. No simulation is run. (Usefull to transport it to another machine)

- **pele_exec**: Use a pele executable that is not the default one. **Needs to be used with pele_data and pele_documents**

- **pele_data**: Use a pele data folder that is not the default one.

- **pele_documents**: Use a pele documents folder that is not the default one.


..  code-block:: yaml

  iterations: 30
  steps: 12
  test: true
  usesrun: false
  debug: true
  pele_exec: "/home/pele/bin/Pele_mpi"
  pele_data: "/home/pele/Data/"
  pele_documents: "/home/pele/Documents/"

Receptor preparation
=======================

Configure the parameters of the PPP (Protein Pele Preparation)

- **preprocess_receptor**: Skip protein pele preparation. Default: False

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


PELE params
================

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



Adaptive params
===================

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

- **metal_constr**: Metal constraints. Default=200

- **water_constr**: Water constraints. Default=5

..  code-block:: yaml

    ca_constr: 2
    interval_constr: 10
    metal_constr: 100
    water_constr: 5


WaterPerturbation
======================

Water modes
+++++++++++++++++

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
++++++++++++++++++++++++

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

- **atom_dist**: Calculate distance between two atomnumbers. Default=None

- **rmsd_pdb**: Calculate rmsd of the ligand to a native pdb structure


..  code-block:: yaml

    atom_dist:
        - 40
        - 1960
    rmsd_pdb: "/home/dsoler/native.pdb"


Output
==========

Configure the output

- **working_folder**: Name of the main working folder where to store the processed input, control files and the simulation folder. Default="resname_Pele_X" where X is a number.

- **output**: Output folder of the simulation. Default=output

..  code-block:: yaml

    working_folder: "NOR_solvent_OBC"
    output: "output_sim"


Automatic Modes
--------------------

Automatically configures all control file options to a standard job chosen beween
induce fit, local exploration, bias exploration, exit path and global exploration


Induced fit
==============

- **induced_fit**: Run induced fit simulation paramaters by setting the center of the box in the
  cm of the ligand, a box radius of 10A, small rotations and translations and a high number of 
  steric clashes and sidechain predition frequency. Usefull to refine docking poses, and search
  new conformations within the same binding site.

..  code-block:: yaml

  induced_fit: true

Rescoring
============

Simulation to refine around an initial conformation. Not looking to find a new binding mode but to minimize
the actual one.

..  code-block:: yaml

  rescoring: true

Local Exploration
=====================

- **out_in**: Local exploration to move the ligand from the bulk to the binding site. The box center set on the 
  center of mass of the ligand with a radius of 30A, steering 1 50% of the times, and a slight bias towards binding energies.
  Useful when no docking is possible in the binding site and you need to open up the pocket.

..  code-block:: yaml

  out_in: true

Biased
=========

- **bias**: Bias exploration towards the indicated bias column. The box center is set on the center of mass of the ligand with
  a radius of 30A, and a bias towards the chosen metric is set. An epsilon fraction of processors are distributed proportionally to the value of a metric, and the rest are inverselyProportional distributed. Therefore, the **epsilon** value controls fraction of the processors that will be assigned according to the selected metric in **biascolumn**


..  code-block:: yaml

  bias: true
  epsilon: 0.5
  bias_column: 5 (starting by 1 on the reports)

Exit path
==============

- **in_out**: Explore the dissociative path of a molecule. At each step the box is center on the most exterior cluster
  and there is a bias towards higher values of SASA. This type accepts a **exit_metric** which represents a column in the report file, an **exit_value** which represents a value for the metric and a **exit_condition** parameter which can be either “<” or “>”, default value is “<”. The simulation will terminate after the metric written in the metricCol reaches a value smaller or greater than exitValue, depending on the condition specified. An example of the exit condition block that would terminate the program after 4 trajectories reaches a value of more than 0.9 for the sixth column (6th starting to count from 1) of the report file would look like:


..  code-block:: yaml

  in_out: true
  exit_value: 0.9
  exit_condition: ">"
  exit_trajnum: 4

Global exploration
=====================

- **global**: Configure a global exploration by randomizing the ligand all around the protein. Then the simulation will start from all configurationsof the system at the same time. The number of configurations (ligand-protein systems) can be chosen thorugh the **poses** flag.

..  code-block:: yaml

  global: true
  poses: 40
