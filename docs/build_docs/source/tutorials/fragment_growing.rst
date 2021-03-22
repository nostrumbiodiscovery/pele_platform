Fragment growing tutorial (unfinished)
========================================

The aim of this tutorial is to familiarise the user with FragPELE, a novel tool to perform fragment growing onto a scaffold
in a high-throughput manner, taking into account the receptor flexibility.

We will be using casein kinase II subunit alpha structure crystallised with a fragment (PDB code 5CSH) as an example system.


1. Protein preparation
+++++++++++++++++++++++

First, we will import the complex from Protein Data Bank and preprocess it using Schrödinger Maestro (release 2020-1).

a. Launch Maestro and **import the structure from PDB** by clicking ``File -> Get PDB...``, type in your PDB ID, e.g. ``5CSH``,
and click ``Download``. The protein structure should appear in your workspace.

.. image:: ../img/fragment_growing_tutorial_1a.png
  :width: 400
  :align: center

b. As you can see, the crystal structure is **duplicated**. Click on ``Select -> Set pick level -> Chains``, then select chain B in the workspace. Open ``Build`` window and remove
selected atoms. This way we will be left with a single protein-scaffold complex.

.. image:: ../img/fragment_growing_tutorial_1b.png
  :width: 400
  :align: center

Similarly, you can select ligand ``54E`` (residue 402) and remove it, since we will not need it.

c. **Preprocess the protein** using Protein Preparation Wizard. Click on ``Tasks`` and search ``Protein Preparation Wizard``.
Check the following options and hit ``Preprocess``.
    - Fill in missing side chain using Prime
    - Fill in missing loops using Prime
    - Delete waters beyond 5.0 Å from het groups

.. image:: ../img/allosteric_tutorial_1b.png
  :width: 400
  :align: center

Maestro will warn you about three residues with double occupancy issues. Simply select them all, click ``Commit`` and then
``Update`` to ensure there are no more problems with the protein.

.. image:: ../img/fragment_growing_tutorial_1c.png
  :width: 400
  :align: center

Upon completion, you should see ``5CSH - preprocessed`` on the entry list.

d. **Change ligand chain ID and residue name**
    - Click on ``Select -> Set pick level -> Residues``, then select ligand ``54E`` (residue 401) with a mouse click
    - Open ``Build`` and choose ``Other edits -> Change atom properties...``
    - Set residue name to ``LIG`` and chain name to ``Z``
    - Choose ``PDB atom name`` from the drop down list and select ``Set unique PDB atom names within residues``
    - Click ``Apply`` and close the window

.. image:: ../img/allosteric_tutorial_1c1.png
  :width: 400
  :align: center

.. image:: ../img/allosteric_tutorial_1c2.png
  :width: 400
  :align: center

e. **Export structure** by clicking on ``File -> Export structures...`` and save all workspace atoms as ``scaffold.pdb``
in your working directory. Close Maestro.

2. Ligands
++++++++++++

FragPELE requires an SD file with fully grown ligands (scaffold + R-group) as input. Maestro provides an easy tool to enumerate
a library of fragments onto a selected hydrogen in your scaffold.

a. **Extract the ligand** to a separate entry by selecting it (just like in step 1d), opening ``Build`` and clicking ``Copy
selected atoms to new entry``.

.. image:: ../img/fragment_growing_tutorial_2a.png
  :width: 400
  :align: center

b. **Enumerate R-groups**

    - Click on ``Select -> Set pick level -> Atoms`` and select any of the hydrogen atoms adjacent to nitrogen
    - Click on ``Tasks``, search for ``Enumeration`` and select ``Enumeration > Custom R-Group...``
    - Choose ``R-groups to Create a Hydrogen Bond`` from the drop down menu and ``Run``

.. image:: ../img/fragment_growing_tutorial_2b.png
  :width: 400
  :align: center

Once the job has finished, you will notice a new group on your entry list called ``enumerate_1``.

c. **Prepare ligands**

Select all enumerated ligands by clicking on the ``enumerate_1`` group, then click on ``Tasks`` and search for ``LigPrep``.
Run LigPrep using the following options (otherwise it will return a lot of structures):

    - use selected structures from the project table
    - generate at most 2 stereoisomers per ligand
    - do not change the ionisation state

.. image:: ../img/fragment_growing_tutorial_2c.png
  :width: 400
  :align: center

d. **Export as SD file**

Once ligand preparation is finished (it might take a few minutes), select all entries in ``ligprep_1-out1`` group and right click on the list
``Export -> Structures...`` and save the file as ``ligands.sdf``.

3. FragPELE input file
++++++++++++++++++++++++

Once you have ``scaffold.pdb`` and ``ligands.sdf`` in your working directory, create an input file called ``input.yaml`` with the following flags:

    - **frag_core** - protein-scaffold PDB file
    - **frag_ligands** - SD file with fully grown ligands
    - **resname** - residue name of the scaffold (we renamed it to ``LIG`` in step 1d)
    - **chain_core** - chain ID of the scaffold (``Z`` in our case)
    - **cpus** - number of CPUs you would like to use for each fragment growing simulation

..  code-block:: yaml

    frag_core: "scaffold.pdb"
    frag_ligands: "ligands.sdf"
    resname: "LIG"
    chain_core: "Z"
    cpus: 50

4. Launching the simulation
+++++++++++++++++++++++++++++

You can launch the simulation from your working directory using one of the following methods:

    - **directly** on command line using ``python -m pele_platform.main input.yaml``

    - submit a slurm file to the **queue system** (ask your IT manager, if you are not sure how to do it). In our case, the slurm file is called ``run.sl`` and we can launch it on the command line using ``sbatch slurm.sl``

Example slurm file:

.. code-block:: console

    #!/bin/bash
    #SBATCH -J PELE
    #SBATCH --output=mpi_%j.out
    #SBATCH --error=mpi_%j.err
    #SBATCH --ntasks=50
    #SBATCH --mem-per-cpu=1000

    python -m pele_platform.main input.yaml

5. Analysis of the results
++++++++++++++++++++++++++++

Once the simulation has finished, your working directory should contain (among other files):

    - a TSV file with scored fragments ``simulation_score.tsv``
    - a number of fragment folders with seemingly random names, e.g. ``scaffold_processed_scaffold,ligand-designer_form-hbond_10evvpbrpvN1-H6C1-H1``

a. Ranked fragments
--------------------

Examine ``simulation_score.tsv`` file

b. Best binding poses
-----------------------

In each fragment folder, you will find a ``top_results`` directory containing PDB files with the best binding poses for that fragment.


