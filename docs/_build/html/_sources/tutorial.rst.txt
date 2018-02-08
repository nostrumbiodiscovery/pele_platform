.. _tutorial:

========
Tutorial
========

1. Setting env variables
----------------------

Before starting, the environment variables the sofware
uses must be set. Therefore, the next commands are needed:

1.1 export set PATH=/path/to/openmpi/binaries/:$PATH

	e.g. $ set export PATH=/usr/lib64/Openmpi/bin:$Path

1.2 export set PATH=/path/to/PELE/binaries:$PATH

	e.g. $ set export PATH=/usr/lib64/Openmpi/bin:$PATH

1.3 export  Schrodinger

	e.g $ set export SCHRODINGER='/opt/schrodinger2016/'


2. Using GUI
-------------

2. python gui.py


2. Using Command Line
----------------------

2. bash [-options] pdbfile ligand_residue ligand_chain
	$ bash PelePlop.sh - -mtor 4  - -clean Samples/jak2_999_complex_processed.pdb LIG Z


