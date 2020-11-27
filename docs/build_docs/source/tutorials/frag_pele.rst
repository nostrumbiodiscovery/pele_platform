FRAG_PELE
========================================

FragPELE is a new tool for in silico hit-to-lead drug design, capable of growing a fragment from a bound core while exploring the protein-ligand conformational space.

This tutorial aims to describe the general protocol to run FragPELE.


Previous Requisites
-----------------------

* **Complex PDB:** The PDB processed file. Prepare the system with the **Schrödinger Protein Prepaation Wizard**. It is recommended to delete water molecules more than 5a away from ligands and ions as well as filling in missii
ng loops and side chains.
It is obligatory that the protein is protonated.
Furthermore, make sure the ligand has:
        * A unique chain ID.
        * Unique PDB atom names with no spaces or single letters.
        * NO residue name except for UNK.

* **Fragment PDB:** The PDB with the desired fragment. The chain of the fragment **must be renamed to "L"**.
* **Serie file:** File with the instructions that explain how the growing is produced and stored. It must have the following format:


Launch a FragPELE simulation
---------------------------------

The default simulation on FragPele is 10 growing steps of 6 Pele steps, and a final simulation of 20 Pele steps.

1. Protein Preparation
-----------------------

	a. Launch Schrodinger Maestro.
 	b. click ``File -> Get PDB`` and type your PDB Id to import the structure.
	c. To preprocess the protein, go to ``Tasks`` and search for ``Protein Preparation Wizard``. Select the following options:
	   Click ``Preprocess`` to start the preprocessing of the protein. 
	d. Change the ligand chain ID ans the residue name.
		#. Go to ``Select -> Set pick level -> Residues``.
		#. Select the ligand with a mouse click.
		#. Go to ``Build`` and click on `` Other edits -> Change atom properties``.
		#. Change ``Residue Name`` to ``LIG``.
		#. Change ``Chain Name`` to ``Z``.
		#. Select ``PDB atom name`` from the drop down list and select ``Set unique PDB atom names within residues``.
		#. Click ``Apply``.
		#. Close the window.
	e. Finally, export the structure by going to ``File -> Export structures`` and save it to your working directory. 

2. Ligand Preparation
------------------------
	a. Select the ligand with a mouse click and extract it to a separate entry opening ``Build`` and clicking ``Copy selected atoms to new entry``. 
	b. Now define the R-groups:
		#. Hit ``Select -> Set pick level -> Atoms``.
		#. Click on nay hydrogen atoms adjacent to Nitrogen.
		#. Go to ``Tasks -> Enumeration -> Custom R-Group``.
		#. Choose ``R-groups to Create a Hydrogen Bond`` from the drop down list.
		#. Click ``Run`` to submit the job. 
	c. An new group on the entry list is created once the job finishes. Select all enumerated ligands by clicking on the group.
	d. Go to ``Tasks -> LigPrep``
	e. Check the following options and hit ``Run``.
	f. A new group on the entry list is created after LigPrep finishes. Select all the netries of the group as in step ``2e``.
	g. Go to ``Export -> Structures`` and save the file as ``ligands.sdf`` in your working directory.

3. YAML Input File
----------------------


4. Launching FragPELE
-----------------------

5. Results
--------------

Optative Flags
------------------

These flags are **exclusive to FragPELE** modes.

- **growing_steps**: Number of steps to grow the fragment with.

- **steps_in_gs**: Number of pele steps within each growing step

- **sampling_steps**: Number of pele steps in the final sampling simulation

- **protocol**: Type of protocol. options = [HT, ES]. For more info please refere here.
        - **HT:** To run FragPELE in **high throughput** mode. 
        - **ES:** 

..  code-block:: yaml

    growing_steps: 6
    steps_in_gs: 6
    sampling_steps: 20
    protocol: HT
    cpus: 24
