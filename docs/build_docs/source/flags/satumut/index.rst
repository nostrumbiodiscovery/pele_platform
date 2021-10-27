Enzyme engineering parameters
=============================

These flags are **common to Saturated mutagenesis and Plurizymer** modes.

- **plurizymer_atom**: Atom to use in the search for neighbouring residues in
  the format chain:resnum:atom_name.
- **satumut_pdb_dir**: The name for the mutated pdb folder (default "pdb_files").
- **satumut_fixed_residues**: List of residues that you don't want to have mutated (Must write the list of residue numbers).
- **satumut_radius_neighbors**: Radius around the *plurizymer_atom* to search
  for residues to mutate (default is 5 angstrom).
- **satumut_hydrogens**: Whether to remove the pdb hydrogen atoms (default False, generraly should be avoided)
- **cpus_per_mutation**: Number of cpus to use for each mutation.

..  code-block:: yaml

    plurizymer_atom: "L:GTP:C4"
    satumut_pdb_dir: "mutated_pdbs"
    satumut_fixed_residues: 
      - 125
      - 180
    satumut_radius_neighbors: 5.3
    satumut_hydrogens: false
    cpus_per_mutation: 5

These flags are **exclusive to saturated mutagenesis** mode. 

- **satumut_positions_mutations**: Residues to mutate, in the format
  chain:resnum (if not specified neighbouring residues around the mutation atom
  will be selected).
- **satumut_mutation**: Residues to mutate to mutate to (using the 3 letter
  code).
- **satumut_multiple_mutations**: Whether to mutate two residues in the same
  pdb (default False).
- **satumut_profile_metric**: The metric to generate the pele profiles with
  (options are *Binding energy* or *currentEnergy*).
- **satumut_consecutive**: Whether to consecutively mutate the PDB file for
  several round (default False).
- **satumut_plots_dpi**: The dpi value to use for the plots.
- **satumut_enantiomer_improve**: The enantiomer that should improve (*R* or *S*).
- **satumut_catalytic_distance**: The distance considered to be catalytic (default is 3.5 angstrom).
- **satumut_dihedrals_analysis**: The four atoms necessary to calculate the
  dihedrals in format chain:resnum:atom_name.
- **satumut_conservative**: How conservative the mutations sould be, options
  are 1 or 2.
- **satumut_energy_threshold**: An energy threshold that limits the points of
  scatter plots.
- **satumut_summary_path**: Name of the summary file created at the end of the
  analysis.
- **satumut_plots_path**: Path of the folder to where to store the plots.
- **satumut_threshold**: The threshold for the improvement which will affect
  what will be included in the summary (default 0). 
- **satumut_analysis_metric**: The metric to measure the improvement of the
  system, options are *energy*, *distance* or *both*.

..  code-block:: yaml

    satumut_positions_mutations: 
      - "A:123"
      - "B:350"
    satumut_mutation: 
      - ALA
      - VAL
      - TYR
    satumut_multiple_mutations: true
    satumut_profile_metric: "currentEnergy"
    satumut_consecutive: false
    satumut_plots_dpi: 600
    satumut_enantiomer_improve: "R"
    satumut_catalytic_distance: 3.5
    satumut_dihedrals_analysis:
      - A:123:H
      - A:123:CA
      - A:123:C
      - A:123:N
    satumut_conservative: 1
    satumut_energy_threshold: 20
    satumut_summary_path: "path/summary"
    satumut_plots_path: "path/plots"
    satumut_threshold: 0.2 
    satumut_analysis_metric: "both"

These flags are **exclusive to plurizymer** mode. 

- **plurizymer_single_mutation**: Name of the residue to mutate to (both 1 and
  3 letter codes can be used).
- **plurizymer_turn**: The number of the round of plurizyme generation, not
  needed for the first round.


..  code-block:: yaml

    plurizymer_single_mutation: "SER"
    plurizymer_turn: 2
    plurizymer_atom: "L:GTP:C4"

