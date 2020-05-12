Versions
############

Here we report the main changes for each version


v1.4.4 (5/05/2020)
=====================

- Include further testing of alignment and rdkit symmetry problem

- Include more flags for FragPele

- Improve exceptions with custom errors

v1.4.3 (27/04/2020)
======================

- Fix rdkit substructure search symmetry problem by alignment

v1.4.2 (23/04/2020)
====================

- FragPELE better tested

- Coverage Platform up to 90%

- Pyyaml checker for unexisting keywords in input

- Improve substructure search on symmetric cases

- Minor fixes

v1.4.1 (23/04/2020)
======================

- Wrongly updated

v1.4.0 (30/03/2020)
=======================

- FragPELE supported (Beta-version)

- PPI simulation supported. Global exploration + induced fit (Beta-version)

- Make Platform work through SCHRODINGER and PELE environment variables

- Get rid of PyMol as external dependency

- Use can define several inputs with asterics. i.e. "complex*.pdb"

- Fix bug on dimer constraints only detecting one chain

- Fix other minor bugs

- Better coverage (77%)


v1.3.4 (10/03/2020)
=======================

- Make mae flag convert clusters as well as top poses to mae

- Let user choose number of clusters through analysis_nclust flag

- Allow user to specify the columns of the report via be_column, te_column and limit_column.

v1.3.3 (01/03/2020)
=======================

- Include only analysis flag

v1.3.2 (30/03/2020)
=======================

- Automatically score the simulation by making the average of the 25% best energy structures.

- Reorder top energy structures

- Support conda deployment for python 3.8

v1.3.1 (29/03/2020)
=======================

- Fixed bug in xtc analysis

- Renew environment on SCHRODINGER subprocess

v1.3.0 (21/02/2020)
=======================

- Set constraints by smiles

- Include a default posprocessing module with plots, top poses and clusters
  
- Separate between AdaptivePELE induced fit (induced_fit_fast) and PELE indeced fit (induced_fit_exhaustive)

- Include skip_ligand_prep option to jump PlopRotTemp missing residue

- Give option ot the user to specify the atom_dist by chain:resname:atomname (A:125:CA)

- Give option mae to transform the best structures to mae files with the metrics as properties

- Fix minor bugs

v1.2.3 (04/02/2020)
=======================

- Automatic PCA mode

- Fix minor bug on global exploration

- Set PPP as external dependence

v1.2.2 (23/12/2019)
=======================

- Fix global exploration bug when joining ligand & receptor

- Add rescoring feature to local a single minimum

- Add induce_fit mode and exploration mode within water_lig parameters to explore hydration sites without moving the ligand or while making the entrance of the ligand.

- Some minor fixes


v1.2.1 (05/12/2019)
=======================

- Add verboseMode

- Add waterPELE and set defaults as we did on WaterMC paper

- Include executable path, data and documents overwriting all constants.py

- Minor fixes

v1.2.0 (24/11/2019)
=======================

- Conda installation

- Insert AdaptivePELE as external dependency

- Fix minor bugs

v1.1.0 (19/10/2019)
=======================

- Automatic Platform to automatically launch PELE&adaptivePELE. It creates the forcefield parameters, the control files, the PELE input.pdb and finally launch the simulation.

- Flexibility to include MSM and Frag PELE

- Flexibility to include analysis scripts

- Flexibility to include PELE modes
