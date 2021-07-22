Understanding the output files
================================

Information about the output of your simulation
-----------------------------------------------------

#. When performing a simulation, a folder will be created. If the simulation is repeated, the folder won't be deleted and a new one will be created.
   The processed input, control files and the simulation folder will be stored there. The default name is: *resname_Pele_X*, where X is a number.

#. The created folder will be organized as:
	* DataLocal (folder): See `PELE Molecular Parameters <https://nostrumbiodiscovery.github.io/pele_docs/molecularParameters.html>`_
	* Lig.mae: Rotamer library
	* docking.pdb: PDB of the docked molecules
	* ligand.pdb: PDB of the ligand
	* receptor.pdb: PDB of the receptor
	* LIG.log: Log file 
	* adaptive.conf: See `Control file outline <https://adaptivepele.github.io/AdaptivePELE/Examples.html#control-file-outline>`_
	* input.yaml: PELE input file
	* pele.conf: Control file. See `General structure of a control file <https://nostrumbiodiscovery.github.io/pele_docs/GeneralStructure/GeneralStructure.html>`_

#. It is recommended to first run the program with the test flag. Thus, the simulation will run using 5 CPUs.

Related Topics
--------------------

* `Installation <../installation/index.html>`_
* `Prepare your own simulation <../packages/index.html>`_
* `Common errors <../errors/index.html>`_
* `Versions <../changelog/index.html>`_
