========================================
FRAG_PELE
========================================

FragPELE is a new tool for in silico hit-to-lead drug design, capable of growing a fragment from a bound core while exploring the protein-ligand conformational space. 

This tutorial aims to describe the general protocol to run FragPELE.

Installation
-----------------------
1. **Conda**
.. code-block:: console
   conda install -c NostrumBiodiscovery -c conda-forge frag_pele --yes

2. **PyPi**

.. code-block:: console
   pip install frag_pele

3. **Source Code**

.. code-block:: console
   git clone https://github.com/carlesperez94/frag_pele.git
   cd frag_pele
   pip install numpy cython #in case not to have them
   python setup.py install


Previous Requisites
-----------------------

* **Complex PDB:** The PDB processed file. Prepare the system with the **Schrödinger Protein Preparation Wizard**. It is recommended to delete water molecules more than 5a away from ligands and ions as well as filling in missing loops and side chains.
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


Flags
-------

1. **High Throughput:** To run in High Throughput use the flap **-HT**. When using this flag FragPELE will perform 3 growing steps of 3 Pele steps anda a final exploration of 20 Pele steps. 

.. code-block:: console
   python frag_pele.main -cp core.pdb -sef serie_file.conf -HT
2. **Standard Precision:** To run in Standard Precision mode us the default values. FragPELE will perform 6 growing steps of 6 Pele steps and a final exploration of 20 Pele steps.

3. **Explorative:** To run in Explorative mode use the flag **-EX**. FragPELE will perform a standard growing simulation with a sampling simulation of:
        * 2 steering
        * 0.5/ 0.3 translation
        * 0.4/0.15 rotation
        * 25 box radius
Select the number of steps in this sampling simulation unsing the flag **-es**.

.. code-block:: console
   python frag_pele.main -cp core.pdb -sef serie_file.conf -sc /path/to/control-personalized.conf
4. **Prepare PDB files:** FragPELE gives the option to only prepare the PDB files. To prepare the PDB files use the flag **-op**.

.. code-block:: console
   python frag_pele.main -cp core.pdb -sef serie_file.conf -op

5. **Only grow fragments:** FragPELE offers the option to only perform the growing of the fragment, if the PDB files were previously prepared. To only execute the growing part of the code use the flag **-og**.

.. code-block:: console
   python frag_pele.main -cp core.pdb -sef serie_file.conf -og

Analysis Tools
------------------

FragPELE also offers several tools to analyse the simulations. All of the analysis tools are available on te path **frag_pele/Analysis**

