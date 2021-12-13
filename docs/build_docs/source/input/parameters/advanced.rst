Advanced parameters
===================

General settings
----------------

Configure the settings of the simulation and the path to all dependencies in case of need (non-default installation).

- **usesrun**: Use srun binary to run PELE. Only when using intel processors.

- **pele_exec**: Use a pele executable that is not the default one. **Needs to be used with pele_data and pele_documents**. default: $PELE/bin/Pele_mpi

- **pele_data**: Use a pele data folder that is not the default one. default: $PELE/Data

- **pele_documents**: Use a pele documents folder that is not the default one. default: $PELE/Documents

- **singularity_exec**: Use a singularity container that contains the pele executable. **Cannot be used with pele_exec**.

- **pele_license**: Use a pele_license path that is not the default one. default: $PELE/licenses

- **schrodinger**: Use a schrodinger path that is not the default one. default: $SCHRODINGER

..  code-block:: yaml

  usesrun: false
  pele_exec: "/home/pele/bin/Pele_mpi"
  pele_data: "/home/pele/Data/"
  pele_documents: "/home/pele/Documents/"
  singularity_exec: "/home/pele/pele_release.sif"
  pele_license: "/home/pele/licenses"
  schrodinger: "/home/pele/schrodinger2020-1/"


PELE general parameters
-----------------------

**These flags are exclusive of the PELE modes not fragPELE**

Advanced parameters of PELE:

- **log**: Retrieve PELE logfiles during simulation. Default=False

- **verbose**: Set to true to activate verbose mode in PELE. DEfault=False

- **anm_displacement**: Angstrom to displace carbon alphas in each ANM movement. Default=0.75

- **anm_modes_change**: Number of steps before we change to a new normal mode movement. Default=4

- **sidechain_res**: Receptor sidechain resolution. Default=10

- **overlap_factor**: Vanderwals overlap factor (More in PELE docs). Default=0.65

- **steric_trials**: Number of steric trials (More in PELE docs). Default=250

- **sidechain_radius**: Residues within specified radius from the ligand will be included in the side chain prediction protocol. Default=6

- **report**: Change the name of the report file. Default: report

- **traj**: Change the name of the trajectory file. Default: trajectory.pdb

..  code-block:: yaml

    log: true
    verbose: true
    anm_displacement: 0.5
    anm_modes_change: 3
    sidechain_res: 30
    overlap_factor: 0.60
    sidechain_radius: 8
    steric_trials: 50
    report: pele_report
    traj: trajectory.xtc


Adaptive PELE parameters
------------------------

Advanced parameters of Adaptive PELE:

- **density**: Density type ([null, exitContinuous...]. More in AdaptivePELE docs). Default: null

..  code-block:: yaml

    density: "exitContinuous"


Preparation parameters
----------------------

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


Residue parametrization
-----------------------

In order to run a simulation, PELE requires the following files for every non-standard molecule (i.e. any non-standard small molecule or residue):

    - **IMPACT template**: containing force field parameters, please refer to `this site <https://nostrumbiodiscovery.github.io/pele_docs/fileFormats.html#impact-template-file-format>`_ for further information.
    - **rotamer library**: optional file containing the list of rotatable bonds to sample by the side chain perturbation algorithm. If missing, the flexibility of the corresponding molecule will not considered. More information available `here <https://nostrumbiodiscovery.github.io/pele_docs/fileFormats.html#sec-fileformats-ligandrotamers>`_.
    - **solvent template**: some special solvents like "OBC" require extra parameters, which are set in this file. Mind that it is **mandatory when using OpenFF**, since it works with the OBC solvent only.

The platform currently has **two implementations** for building hetero molecule parameters - PlopRotTemp (soon to be deprecated) and
`Peleffy <https://github.com/martimunicoy/peleffy>`_ (PELE Force Field Yielder), which offers more functionality but is still in beta testing.

Please refer to the following table for the comparison of the two methods and available forcefields:

+-------------+----------------------+--------------+------------------------------------+
| **Builder** | **Forcefields**      | **Solvents** | **Charge parametrization methods** |
+-------------+----------------------+--------------+------------------------------------+
| PlopRotTemp | "OPLS2005"           | "OBC"        | "OPLS2005"                         |
|             |                      |              |                                    |
| (default)   |                      | "VDGBNP"     |                                    |
+-------------+----------------------+--------------+------------------------------------+
| Peleffy     | "OPLS2005" (default) | "OBC"        | "gasteiger"                        |
|             |                      |              |                                    |
| (beta)      | "openff-1.3.0"       | "VDGBNP"     | "am1bcc" (default for OpenFF)      |
|             |                      |              |                                    |
|             | "openff-1.2.1"       |              | "OPLS2005" (default for OPLS2005)  |
|             |                      |              |                                    |
|             | "openff-1.2.0"       |              |                                    |
|             |                      |              |                                    |
|             | "openff-1.1.1"       |              |                                    |
|             |                      |              |                                    |
|             | "openff-1.1.0"       |              |                                    |
|             |                      |              |                                    |
|             | "openff-1.0.1"       |              |                                    |
|             |                      |              |                                    |
|             | "openff-1.0.0"       |              |                                    |
+-------------+----------------------+--------------+------------------------------------+


PlopRotTemp
+++++++++++

To continue using PlopRotTemp, you do not need to make any changes to your YAML file, previously existing flags are still
available:

    - **gridres**: Resolution of the rotamers when sampling them by the Side Chain prediction algorithm. Default=10 degrees

    - **core**: List of PDB atom names that will be included as part of the rigid core. In case it is not specified, the algorithm will pick up a set of non-rotatable atoms centered in the molecular structure. Default=None

    - **exclude_terminal_rotamers**: Exclude terminal rotamers during parametrization of non standard molecules if they belong to a small terminal group. Default=True

    - **mae_lig**: External MAE file with quantum charges generated with Schrödinger suite. When supplied, any charge calculated internally in the platform will be replaced by the charges from this file. Default=None

    - **maxtorsion**: Maximum number of rotamers per flexible side chain. Default=4

    - **n**: Maximum number of flexible side chains in a molecule. Default=None


..  code-block:: yaml

    solvent: "OBC"
    maxtorsion: 4
    n: 5
    mae_lig: "/home/dsoler/lig.mae"
    gridres: 10


Peleffy
+++++++

In order to use Peleffy instead of PlopRotTemp, you need to set ``use_peleffy: true`` in input YAML.

You can use the following parameters to control the way peleffy will parametrize non-standard molecules for you:

- **forcefield**: Forcefield used to parametrize hetero molecules, you can use one of:

        - "OPLS2005" (default)
        - "openff-1.3.0"
        - "openff-1.2.1"
        - "openff-1.2.0"
        - "openff-1.1.1"
        - "openff-1.1.0"
        - "openff-1.0.1"
        - "openff-1.0.0"

- **charge_parametrization_method**: The method to use to assign partial charges to atoms:

        - "gasteiger"
        - "am1bcc" (default when using any "OpenFF" force field)
        - "OPLS2005" (default when using "OPLS2005")

- **use_peleffy**: You have to set it to True to use peleffy instead of the default parameters builder. Default=False

- **gridres**: Resolution of the rotamers when sampling them by the Side Chain prediction algorithm. Default=10 degrees

- **core**: List of PDB atom names that will be included as part of the rigid core. In case it is not specified, the algorithm will pick up a set of non-rotatable atoms centered in the molecular structure. Default=None

- **exclude_terminal_rotamers**: Exclude terminal rotamers during parametrization of non standard molecules if they belong to a small terminal group. Default=True

- **mae_lig**: External MAE file with quantum charges generated with Schrödinger suite. When supplied, any charge calculated internally in the platform will be replaced by the charges from this file. Default=None

Important: Peleffy requires CONECT lines in the PDB file, otherwise they are automatically added with Schrödinger Protein Preparation Wizard.

..  code-block:: yaml

    use_peleffy: true
    charge_parametrization_method: "gasteiger"
    forcefield: "openff-1.3.0"
    gridres: 20
    core:
        - "O1"
        - "C1"
        - "C2"
        - "N1"


Ligand parameters
-----------------

Advanced parameters for the conformation perturbation algorithm:

- **overlap_factor_conformation**: van der Waals overlap factor in conformation perturbation. Default = 0.65

..  code-block:: yaml

    overlap_factor_conformation: 0.60


Constraint parameters
+++++++++++++++++++++

Alternatively to constraint levels, advanced users can manipulate the constraint parameters individually at their own risk, using the following flags:

- **terminal_constr** - sets the spring constant for the terminal C-alpha constraints, default = 5 kcal/mol

- **ca_constr** - sets the spring constant for the remaining C-alphas in the backbone, default = 0.5 kcal/mol

- **ca_interval** - interval at which the backbone C-alphas should be constrained, default = 10 (i.e. every 10 residues).

Take into account that specific modifiers of constraint parameters will prevail over the settings coming from the
constraints levels and those predefined in each package.

Other advanced parameters related with constraints:

- **water_constr**: Water constraints. Default=5

- **constrain_smiles**: SMILES string to indicate what part of the molecule to constrain. Default=None

- **external_constraints**: You can specify 2 types of constraints: positional constraints or atom-atom constraints, e.g.

  - The positional constraints are given either by:
        - springConstant-atomnumber. i.e. "10-17"
        - springConstant-chain:resnum:atomname. i.e. "5-A:1:H"

  - The atom-atom constraints are specified either by:
        - springConstant-equilibriumDistance-atomnumber1-atomnumber2. i.e. "50-2.34-17-4159"
        - springConstant-equilibriumDistance-chain1:resnum1:atomname1-chain2:resnum2:atomname2. i.e. "50-2.34-A:1:H-L:1:C21"

- **remove_constraints**: Do not place constraints on the carbon-alpha of the protein. Default: False

..  code-block:: yaml

    terminal_constr: 10.5
    ca_constr: 6.0
    ca_interval: 3
    water_constr: 5
    constrain_smiles: "C2CCC1CCCCC1C2"
    smiles_constr: 5
    external_constraints:
    - "10-17" #constrain of 10kcal/mol at atomnumber 17
    - "5-A:1:H" ##constrain of 10kcal/mol at atom with chain A residuenumber 1 and atomname H
    - "50-2.34-17-4159" #constrain of 50kcal/mol with equilibrium distance of 2.34 between atomnumbers 17 & 4159
    - "50-2.34-A:1:H-L:1:C21" #constrain of 50kcal/mol with equilibrium distance of 2.34 between atoms with respective chain resnum and atomname
    remove_constraints: true


Metal constraints
+++++++++++++++++

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
    constrain_core: "CN(C)C(=O)c1ccc(F)cc1"  # SMILES or SMARTS pattern
    constrain_core_spring: 30  # optional, default 50.0


aquaPELE parameters
-------------------

Other advanced parameters to set up aquaPELE. Its usage is discouraged.

- **water_trials**: Numerical trials on water perturbation. Default=10000

- **water_constr**: COM constrain applied to th water molecule after perturbation. Default=0

- **water_overlap**: Overlap factor of water. Default=0.78


..  code-block:: yaml

    water_trials: 500
    water_constr: 0.5
    water_overlap: 0.5


Analysis parameters
-------------------

Advanced parameters of Analysis package.

- **analysis_nclust**: Numbers of clusters out of the simulation, if using the ``gaussianmixture`` clustering method. Default: 10

- **be_column**: Column of the binding energy in the reports starting by 1. Default: 5

- **te_column**: Column of the total energy in the reports starting by 1. Default: 4

- **limit_column**: Specify the column where your external metrics start. Default: 6

- **mae**: To extract the best energy and cluster poses as .mae files with the metrics as properties (schrodinger need it). Default: false

- **clustering_method**: If you want to override the default clustering method (meanshift), you can set this flag to ``gaussianmixture`` or ``HDBSCAN``.

- **clustering_filtering_threshold**: Percentage of output structures to filter our before clustering. Default = 0.25.

- **plot_filtering_threshold**: Percentage of output structures to filter out before creating plots. Default = 0.02

..  code-block:: yaml

    analysis_nclust: 12
    be_column: 5
    te_column: 4
    limit_column: 6
    mae: true
    clustering_method: "gaussianmixture"
    clustering_filtering_threshold: 0.1
    plot_filtering_threshold: 0.1

.. note::
   In case of the ``meanshift`` algorithm,
   the `bandwidth <basic_parameters/analysis.html#bandwidth>`__ refers to the
   maximum RMSD allowed within the cluster, whereas in ``HDBSCAN`` to distances
   between your data points.


Output
------

Configure the output

- **output**: Output folder of the simulation. Default=output

..  code-block:: yaml

    output: "output_sim"
