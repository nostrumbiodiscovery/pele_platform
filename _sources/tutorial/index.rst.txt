Launch your first job
##########################

0) Prepare the system with maestro (Protein Preparation Wizard, hydrogen optimization and posterior minimization)
and output a complex.pdb

1) Copy the text below into a file called input.yaml and change the path
to your input pdb file, the residue name of the ligand to be perturb and its chain.

..  code-block:: yaml

    system: "complex.pdb"
    residue: "LIG"
    chain: "L"

2)  Choose one of the modes below by copying its line to the input.yaml:

..  code-block:: yaml

    induced_fit: true
    out_in: true
    in_out: true
    global: true   
    
3) Test your simulation (1min):

First add the flag **test:true** into and launch the next command:

``python -m pele_platform.main input.yaml``

4) Run prodution

If everything works change the flag **test: true** by **cpus: X** being X the number of cpus to use for the job. From 40-60 for induced fit simulation and exit paths, from 100-200 for local explorations and 250 for global. Finally run the python command again:

``python -m pele_platform.main input.yaml``

