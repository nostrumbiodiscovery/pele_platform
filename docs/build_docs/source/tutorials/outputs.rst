*****************
Output files
*****************




Considerations when using Pele Platform
###########################################

1. When performing a simulation, a folder will be created. If the simulation is repeated, the folder won't be deleted and a new one will be created.
   The processed input, control files and the simulation folder will be stored there. The default name is: *resname_Pele_X*, where X is a number.
2. The created folder will be organized as:
	* DataLocal (folder): See `PELE Molecular Parameters <https://eapm-bsc.github.io/PELE-repo/molecularParameters.html>`_ 
	* Lig.mae: Rotamer library
	* docking.pdb: PDB of the docked molecules
	* ligand.pdb: PDB of the ligand
	* receptor.pdb: PDB of the receptor
	* LIG.log: Log file 
	* adaptive.conf: See `Control file outline <https://adaptivepele.github.io/AdaptivePELE/Examples.html#control-file-outline>`_
	* input.yaml: PELE input file
	* pele.conf: Control file. See `General structure of a control file <https://eapm-bsc.github.io/PELE-repo/GeneralStructure/GeneralStructure.html>`_
3. It is recommended to first run the program with the test flag. Thus, the simulation will run using 5 CPUs.

Related Topics
#################

* `Installation <../installation/index.rst>`_
* `Prepare your own simulation <../packages/index.rst>`_
* `Common errors <../errors/index.rst>`_
* `Versions <../changelog/index.rst>`_

.. toctree::
   installation/index.rst
   :hidden::
.. toctree::
   packages/index.rst
   :hidden::
.. toctree::
   errors/index.rst
   :hidden::
.. toctree::
   changelog/index.rst
   :hidden::


