Versions
############

Here we report the main changes for each version


v1.3
=======================

- Set constraints by smiles

- Include a default posprocessing module with plots, top poses and clusters
  
- Separate between AdaptivePELE induced fit (induced_fit_fast) and PELE indeced fit (induced_fit_exhaustive)

- Include skip_ligand_prep option to jump PlopRotTemp missing residue

- Give option ot the user to specify the atom_dist by chain:resname:atomname (A:125:CA)

- Give option mae to transform the best structures to mae files with the metrics as properties

- Fix minor bugs

v1.2.3
=======================

- Automatic PCA mode

- Fix minor bug on global exploration

- Set PPP as external dependence

v1.2.2
=======================

- Fix global exploration bug when joining ligand & receptor

- Add rescoring feature to local a single minimum

- Add induce_fit mode and exploration mode within water_lig parameters to explore hydration sites without moving the ligand or while making the entrance of the ligand.

- Some minor fixes


v1.2.1
=======================

- Add verboseMode

- Add waterPELE and set defaults as we did on WaterMC paper

- Include executable path, data and documents overwriting all constants.py

- Minor fixes

v1.2.0
=======================

- Conda installation

- Insert AdaptivePELE as external dependency

- Fix minor bugs

v1.1
=======================

- Automatic Platform to automatically launch PELE&adaptivePELE. It creates the forcefield parameters, the control files, the PELE input.pdb and finally launch the simulation.

- Flexibility to include MSM and Frag PELE

- Flexibility to include analysis scripts

- Flexibility to include PELE modes
