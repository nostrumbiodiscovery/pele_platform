Fragment growing simulation
===============================

Introduction
----------------

FragPELE is a new tool for in silico hit-to-lead drug design, capable of performing high-throughput fragment growing
onto a scaffold while exploring the protein-ligand conformational space. It aims to rank ligands and their binding modes,
as well as discover cryptic pockets.

For more details, you can refer to:

- `FragPELE: Dynamic Ligand Growing within a Binding Site. A Novel Tool for Hit-To-Lead Drug Design <https://pubmed.ncbi.nlm.nih.gov/32027130/>`_
- `frag_pele documantation <https://carlesperez94.github.io/frag_pele/>`_.

Inputs
++++++++++++++
- protein-scaffold PDB file
- grown ligands SD file
- input.yaml configuration file

Default parameters
+++++++++++++++++++

- frag_eq_steps: 20
- gr_steps: 6
- frag_steps: 3

Recommendations
+++++++++++++++++

#. We recommend using **at least 50 CPUs**.
#. The simulation can take between 30 mins and 3 hours per fragment, depending on the number of rotatable bonds.

1. Complex Preparation
------------------------
   
Protein-scaffold PDB file
++++++++++++++++++++++++++++++

PDB file should be preprocessed with Maestro Protein Preparation Wizard. We usually recommend protonating the protein (obligatory), deleting water molecules more than 5Å away from ligands
and ions as well as filling in missing loops and side chains.

The file must contain protein in complex with a fragment (scaffold), e.g. an X-ray structure or a docked pose. Additionally, ensure the ligand has:

 - unique chain ID
 - no atom names with spaces or single letters (occasionally Maestro adds hydrogens named ``H1 2``, these need to be corrected)
 - any unique residue name, except for ``UNK``
 - well-defined aromatic bonds.

SD file with fully grown ligands
++++++++++++++++++++++++++++++++++

Ligands should be preprocessed using Schrodinger LigPrep (default settings should be sufficient) and have unique molecule names.

2. Input Preparation
-----------------------
 
Prepare the input file ``input.yml``:

    - **frag_core** - path to preprocessed PDB file containing the protein and docked scaffold
    - **frag_ligands** - path to SD file with fully grown, preprocessed ligands
    - **resname** - unique residue name of the scaffold
    - **chain_core** - unique chain ID of the scaffold
    - **cpus** - number of CPUs to use

..  code-block:: yaml

    frag_core: "/home/daniel/PR_core.pdb" # Protein-scaffold complex
    frag_ligands: "/home/daniel/grown_ligands.sdf" # Fully grown ligands
    resname: "LIG" # Ligand scaffold residue name
    chain_core: "L" # Ligand scaffold chain ID
    cpus: 48

For more optional flags please refer to `optional flags <../../flags/index.html>`_.


3. Run simulation
---------------------

To run the system launch the simulation with the following command:

``python -m pele_platform.main input.yml``

4. Output
--------------

The simulation will create a TSV file with scored fragments as well as a number of fragment folders, the names of which will consist of the scaffold and molecule names.

Scored fragments
+++++++++++++++++++++

The list of all grown fragments together with their associated scores (average binding energy of the top 25% of all poses)
can be found in:

``simulation_score.tsv``


Top poses
+++++++++++++

Each fragment folder contains a top_results folder with PDB files corresponding to the best poses for that fragment:

``fragment_folder/top_results/``
