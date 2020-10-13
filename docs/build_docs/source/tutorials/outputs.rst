*****************
Output files
*****************




Considerations when using Pele Platform
###########################################

#. When performing a simulation, a folder will be created. If the simulation is repeated, the folder won't be deleted and a new one will be created.
   The processed input, control files and the simulation folder will be stored there. The default name is: *resname_Pele_X*, where X is a number.

#. The created folder will be organized as:
	* DataLocal (folder): See `PELE Molecular Parameters <https://eapm-bsc.github.io/PELE-repo/molecularParameters.html>`_ 
	* Lig.mae: Rotamer library
	* docking.pdb: PDB of the docked molecules
	* ligand.pdb: PDB of the ligand
	* receptor.pdb: PDB of the receptor
	* LIG.log: log file 
	* adaptive.conf: Control file outline
	* input.yaml: PELE input file
	* pele.conf: Control file
#. It is recommended to first run the program with the test flag. Thus, the simulation will run using 5 CPU's.

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


