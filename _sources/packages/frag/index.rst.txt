Prepare your own Fragment Growing
######################################

FragPELE is a new tool for in silico hit-to-lead drug design, capable of performing HT-fragment growing onto a scaffold while
exploring the protein-ligand conformational space. Aims to rank ligands and discover cryptic pockets.

**Article**: https://nostrumbiodiscovery.github.io/papers/Software/index.html#frag-pele-dynamic-ligand-growing-within-a-binding-site-a-novel-tool-for-hit-to-lead-drug-design

**Input** (further explained below):


    - protein-scaffold PDB file
    - grown ligands SD file
    - input.yaml configuration file

**Output**:

    - ranking of the fragments
    - ranked binding modes within each fragment

**Computational time**:

    - 30min/3h per fragment (depends on how many rotatable bonds the fragment has)

1. Complex Preparation
======================
   
**Protein-scaffold PDB file** should be preprocessed with Maestro Protein Preparation Wizard.
We would usually recommend protonating the protein (obligatory), deleting water molecules more than 5Å away from ligands
and ions as well as filling in missing loops and side chains.

The file must contain protein in complex with a fragment (scaffold), e.g. an X-ray structure or a docked pose. Additionally, ensure the ligand has:

 - unique chain ID
 - no atom names with spaces or single letters (occasionally Maestro adds hydrogens named ``H1 2``, these need to be corrected)
 - any unique residue name, except for UNK
 - well-defined aromatic bonds.

**SD file with fully grown ligands** which should be preprocessed using Schrodinger LigPrep (default settings should be sufficient) and have unique molecule names.

2. Input Preparation
=====================
 
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

For more optional flags please refer to `optative falgs <../../documentation/index.html>`_


3. Run simulation
====================


To run the system launch the simulation with the following command:

``python -m pele_platform.main input.yml``

4. Output
===============

The simulation will create a TSV file with scored fragments as well as a number of fragment folders, the names of which will consist of the scaffold and molecule names.

**Scored fragments**

The list of all grown fragments together with their associated scores (average binding energy of the top 25% of all poses)
can be found in:

``simulation_score.tsv``


**Top poses**

Each fragment folder contains a top_results folder with PDB files corresponding to the best poses for that fragment:

``fragment_folder/top_results/``


