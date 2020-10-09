*****
Output files
*****




Considerations when using Pele Platform
#####

#. When performing a simulation, a folder will be created. If the simulation is repeated, the folder won't be deleted and a new one will be created.
   The processed input, ontrol files and the simulation folder will be stored there. The default name is: *resname_Pele_X*, where X is a number.

#. The created folder will be organizes as:
	* DataLocal (folder): See PELE Molecular Parameters
	* Lig.mae: Rotamer library
	* docking.pdb: PDB of the docked molecules
	* ligand.pdb: PDB of the ligand
	* receptor.pdb: PDB of the receptor
	* LIG.log: 
	* adaptive.conf: Control file outline
	* input.yaml
	* pele.conf
#. It is recommended to first run the program with the test flag. Thus, the simulation will be run using 5 CPU's.

Related Topics
#####

* `Installation <installation/index.rst>`
* Prepare your own simulation
* Common errors
* Versions
