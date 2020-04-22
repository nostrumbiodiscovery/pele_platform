Prepare your own PELE simulation
####################################

To prepare your own montecarlo simulation to explore the  protein-ligand conformational space
follow the steps below.

1. Complex Preparation
======================
   
Prepare the system with maestro (Protein Preparation Wizard)
and output a complex.pdb. The complex.pdb must contain the protein-ligand in the desired initial configuration.
**If using the pocket exploration mode, the ligand can be place anywhere**

Make sure the ligand has:

 - Unique chain
 - No atomnames with spaces or single letter
 - Any residuename except UNK

2. Input Preparation
=====================
 
Prepare the input file ``input.yml``:

To run different modes prepare different control files.
For more explanation on the modes please refer to `here <../../modes/adaptive/index.html>`__


Induced fit fast (3-4 h)
+++++++++++++++++++++++++++++++++++++++++

..  code-block:: yaml

    system: 'docking2grid6n4b_thc.pdb'
    chain: 'L'
    resname: 'THC'
    seed: 12345
    induced_fit_fast: true
    atom_dist:
    - "A:2:CA" #First atom to make the distance to
    - "B:3:CG" #Second atom to make the distance to 
    cpus: 60

Induced fit exhaustive (12h)
+++++++++++++++++++++++++++++++++++++++++

..  code-block:: yaml

    system: 'docking2grid6n4b_thc.pdb'
    chain: 'L'
    resname: 'THC'
    seed: 12345
    induced_fit_exhaustive: true
    atom_dist:
    - "A:2:CA" #First atom to make the distance to
    - "B:3:CG" #Second atom to make the distance to 
    cpus: 60

Pose minimization (30min)
+++++++++++++++++++++++++++

..  code-block:: yaml

    system: 'docking2grid6n4b_thc.pdb'
    chain: 'L'
    resname: 'THC'
    seed: 12345
    rescoring: true
    atom_dist:
    - "A:2:CA" #First atom to make the distance to
    - "B:3:CG" #Second atom to make the distance to 
    cpus: 60


Pocket Exploration (24h)
+++++++++++++++++++++++++++

..  code-block:: yaml

   system: 'docking2grid6n4b_thc.pdb'
   chain: 'L'
   resname: 'THC'
   seed: 12345
   atom_dist:
   - "A:2:CA" #First atom to make the distance to
   - "B:3:CG" #Second atom to make the distance to 
   ppi: true

Water+ligand exploration (4h)
++++++++++++++++++++++++++++++++++++

..  code-block:: yaml

    system: 'docking2grid6n4b_thc.pdb'
    chain: 'L'
    resname: 'THC'
    seed: 12345
    water_lig:
    - M:1 #Chain and residue of 1st water
    - M:2 #Chain and residue of 2nd water
    atom_dist:
    - "A:2:CA" #First atom to make the distance to
    - "B:3:CG" #Second atom to make the distance to 
    cpus: 128
    
Out_in local exploration (12h)
++++++++++++++++++++++++++++++++++++

..  code-block:: yaml

    system: 'docking2grid6n4b_thc.pdb'
    chain: 'L'
    resname: 'THC'
    seed: 12345
    out_in: true
    atom_dist:
    - "A:2:CA" #First atom to make the distance to
    - "B:3:CG" #Second atom to make the distance to 
    cpus: 200

Biased exploration (6h)
++++++++++++++++++++++++++

..  code-block:: yaml

    system: 'docking2grid6n4b_thc.pdb'
    chain: 'L'
    resname: 'THC'
    seed: 12345
    bias: true
    epsilon: 0.5
    bias_column: 5 #starting by 1 on the reports
    atom_dist:
    - "A:2:CA" #First atom to make the distance to
    - "B:3:CG" #Second atom to make the distance to 
    cpus: 200

Receptor sampling simulation (10 min)
++++++++++++++++++++++++++++++++++++++++

..  code-block:: yaml

   system: 'docking2grid6n4b_thc.pdb'
   chain: 'L'
   resname: 'THC'
   seed: 12345
   spawning: independent
   ca_constr: 3
   pca_traj:
   - "pele_platform/Examples/pca/1.pdb"
   - "pele_platform/Examples/pca/2.pdb"
   - "pele_platform/Examples/pca/3.pdb"
   selection_to_perturb: false
   perturbation: false
   remove_constraints: true
   cpus: 20
   iterations: 1
   steps: 20

3. Run simulation
====================


To run the system launch the simulation with the next command:

``python -m pele_platform.main input.yml``

