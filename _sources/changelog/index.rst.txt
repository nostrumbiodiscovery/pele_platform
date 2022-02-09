===============
Release history
===============


PELE Platform 1.6.3
-------------------

A minor stability release of the PELE Platform that contains small bug fixes for different packages and compatibility changes to support future PELE versions. Find below the full list of modifications:

- Minor changes to support new beta version of PELE (1.7.2).
- Fixed bug in equilibration_steps parameter.
- Fixed bug in frag PELE.
- Minor documentation modifications.


PELE Platform 1.6.2
-------------------

A minor release of the PELE Platform that introduces small improvements and more flexibility to handle different MPI architectures. Find below the full list of modifications:

- Support for custom MPI parameters.
- New ``sidechain_radius`` flag to control sidechain prediction region radius.
- Meaningful names to describe atom_dist metrics.
- Upgrade peleffy to 1.4.1


PELE Platform 1.6.1
-------------------

A minor release of the PELE Platform that contains stability and efficiency improvements. Find below the full list of modifications:

- Updated randomization of the initial ligand poses, so that they will only be spawned in a user-defined box. As a result, the pocket search can be focused on a specific part of the protein.

- Introduced automated bandwidth selection for the meanshift algorithm.

- PELE Plotter 1.1 with new plot types and filtering option

- Saving analysis parameters to file for later inspection.

- Added support for Open Force Field 2.0.

- Introduced automated selection of cluster conditions for AdaptivePELE with the new preequilibration stage.

- Induced fit exhaustive is renamed to induced fit long.

- Fixed SameFileError exception raised during ligand parametrization.

- Raised an exception when ligand parametrization fails (previously a warning), so it does not fail silently.

- Fixed missing logger bugs when restarting SiteFinder and SaturatedMutagenesis.

- Froze pyparsing version to avoid dependency clashes with prody.

- Automated removal of capped termini.

- Changed default parameters for the OutIn package.

- Updates to documentation, including OutIn and ligand parametrization tutorials.


PELE Platform 1.6
-----------------

The new release of the PELE Platform contains many new features and stability changes. The most important ones are:

- **Fragment libraries**: provide us with a PDB or SD file containing your fragments library and we will enumerate all possible ligands, grow them inside your binding site and rank them from best to worst, always taking receptor flexibility into account.

- **Conformation perturbation**: if you provide PELE with a conformation library, it will perform an additional perturbation step focused on changing the conformation. Especially useful if you want to focus on bioactive conformers only or to apply specific sampling to tricky compounds like macrocycles or heterocyclic compounds.

- **Brand new analysis**: we give users more flexibility and tools to simplify the analysis of PELE simulations. You can now choose your preferred clustering method, its parameters and criteria to extract most promising structures. We also improved the plots and added information about the water clusters resulting from aquaPELE simulations.

- **Covalent docking**: we have developed a covalent docking protocol to let you explore all those irreversibly binding ligands. All you have to do is provide the platform with a ligand covalently attached to the required protein residue and it will sample local conformations to optimize the interactions for you.

- **New built-in Plotter**: a new plotter to quickly represent the results of your simulation using the command line. You can choose from a number of customizable scatter, density and interactive plots.

- **Interaction restrictions**: with the new release, you will be able to restrict the interactions between the ligand and the protein by specifying a custom H-bond distance or angle range between selected atoms.

- **Open Force Field integration**: from now on you can parametrize your ligands and cofactors using Open Force Field as well as the good old OPLS2005.

- **Singularity**: from now on, we support launching PELE from Singularity containers to save some of the trouble with dependencies and compatibility.

Other changes are:

- Kernel density estimator plots

- Constraints levels for alpha carbons

- Changes to Monte Carlo parameters

- Allosteric package renamed to site_finder

- Minor changes to folder structure

- Interstep logger integration

- New optional flags to control top clusters and their representative structures selection

- Recover restart flag to allow the users to manually curate control files and use them in a new run

- Support for the new equilibration mode flag of Adaptive

- Fix problems with the box of the Out --> In package

- Add user warnings to facilitate the system preparation

- Add water sites analysis

- New environment variables PELE_EXEC, PELE_DATA and PELE_DOCUMENTS

- Fix for the site_finder randomization

- Support for a new PELE flag called minimum steps

- Improvements in the docs

- Version header

- New checker for flag incompatibilities

- New checker for input PDB files

- Multi representative structures option

- Frag-3.1.1 upgrade

- Fix for atom_dist bug in biased simulations


v1.5.1
------

- AquaPELE

- High-throughput fragment screening

- Improved out-in exploration

- Support for non-standard residues

- Automatic metal constraints

- Metal polarisation

- Tutorials

- Outliers removed from plots

- Improved documentation


v1.5.0
------

- PPI package

- Site finder (pocket exploration) package

- GPCR orthosteric package

- Binding package

- External metal constraints

- Add n random water to your simulation by setting n_waters flag

- More robust error handling

- Remove support python 3.6 and update features python 3.7

- Full refactor of code

- Improvement of frag_pele

- New docs

- Coverage up to 94%


v1.4.4
------

- Include further testing of alignment and rdkit symmetry problem

- Include more flags for FragPele

- Improve exceptions with custom errors


v1.4.3
------

- Fix rdkit substructure search symmetry problem by alignment


v1.4.2
------

- FragPELE better tested

- Coverage Platform up to 90%

- Pyyaml checker for unexisting keywords in input

- Improve substructure search on symmetric cases

- Minor fixes


v1.4.1
------

- Wrongly updated


v1.4.0
------

- FragPELE supported (Beta-version)

- PPI simulation supported. Global exploration + induced fit (Beta-version)

- Make Platform work through SCHRODINGER and PELE environment variables

- Get rid of PyMol as external dependency

- Use can define several inputs with asterics. i.e. "complex*.pdb"

- Fix bug on dimer constraints only detecting one chain

- Fix other minor bugs

- Better coverage (77%)


v1.3.4
------

- Make mae flag convert clusters as well as top poses to mae

- Let user choose number of clusters through analysis_nclust flag

- Allow user to specify the columns of the report via be_column, te_column and limit_column.


v1.3.3
------

- Include only analysis flag


v1.3.2
------

- Automatically score the simulation by making the average of the 25% best energy structures.

- Reorder top energy structures

- Support conda deployment for python 3.8


v1.3.1
------

- Fixed bug in xtc analysis

- Renew environment on SCHRODINGER subprocess


v1.3.0 
------

- Set constraints by smiles

- Include a default posprocessing module with plots, top poses and clusters
  
- Separate between AdaptivePELE induced fit (induced_fit_fast) and PELE indeced fit (induced_fit_long)

- Include skip_ligand_prep option to jump PlopRotTemp missing residue

- Give option ot the user to specify the atom_dist by chain:resname:atomname (A:125:CA)

- Give option mae to transform the best structures to mae files with the metrics as properties

- Fix minor bugs


v1.2.3
------

- Automatic PCA mode

- Fix minor bug on global exploration

- Set PPP as external dependence


v1.2.2
------

- Fix global exploration bug when joining ligand & receptor

- Add rescoring feature to local a single minimum

- Add induce_fit mode and exploration mode within water_lig parameters to explore hydration sites without moving the ligand or while making the entrance of the ligand.

- Some minor fixes


v1.2.1
------

- Add verboseMode

- Add waterPELE and set defaults as we did on WaterMC paper

- Include executable path, data and documents overwriting all constants.py

- Minor fixes


v1.2.0
------

- Conda installation

- Insert AdaptivePELE as external dependency

- Fix minor bugs


v1.1.0
------

- Automatic Platform to automatically launch PELE&adaptivePELE. It creates the forcefield parameters, the control files, the PELE input.pdb and finally launch the simulation.

- Flexibility to include MSM and Frag PELE

- Flexibility to include analysis scripts

- Flexibility to include PELE modes
