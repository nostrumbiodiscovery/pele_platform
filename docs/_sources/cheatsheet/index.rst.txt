Prepare your own simulation
####################################

0. Prepare the system with maestro (Protein Preparation Wizard, hydrogen optimization and posterior minimization)
and output a complex.pdb.

Make sure the ligand has:

 - Unique chain
 - No atomnames with spaces or single letter
 - Any residuename except UNK

1. Prepare the input file ``input.yml``:

To run different modes prepare different control files

a) **Induce fit**:

..  code-block:: yaml

    system: 'docking2grid6n4b_thc.pdb'
    chain: 'L'
    resname: 'THC'
    seed: 12345
    induced_fit: true
    atom_dist:
    - 4059 #atomnumber1
    - 4556 #atomnumber2
    spawning: independent #usePELEnotAdaptivePELE
    cpus: 200
    iterations: 1
    steps: 1000

b) **Global exploration**:

..  code-block:: yaml

    system: 'docking2grid6n4b_thc.pdb'
    chain: 'L'
    resname: 'THC'
    seed: 12345
    global: true
    poses: 120
    atom_dist:
    - 4059
    - 4556
    cpus: 250

c) **Water+ligand exploration**:

..  code-block:: yaml

    system: 'docking2grid6n4b_thc.pdb'
    chain: 'L'
    resname: 'THC'
    seed: 12345
    water_lig:
    - M:1 #Chain and residue of 1st water
    - M:2 #Chain and residue of 2nd water
    atom_dist:
    - 4059
    - 4556
    cpus: 128
    
d) **Out_in local exploration**:

..  code-block:: yaml

    system: 'docking2grid6n4b_thc.pdb'
    chain: 'L'
    resname: 'THC'
    seed: 12345
    out_in: true
    atom_dist:
    - 4059
    - 4556
    cpus: 200

e) **Biased exploration**:

..  code-block:: yaml

    system: 'docking2grid6n4b_thc.pdb'
    chain: 'L'
    resname: 'THC'
    seed: 12345
    bias: true
    epsilon: 0.5
    bias_column: 5 #starting by 1 on the reports
    atom_dist:
    - 4059
    - 4556
    cpus: 200


2. Run simulation

``python -m pele_platform.main input.yml``

