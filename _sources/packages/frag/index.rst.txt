Prepare your own Fragment Growing
######################################

FrAG, a new tool for in silico hit-to-lead drug design, capable of performing HT-fragment growing into a core while exploring the protein-ligand conformational space. Aims to rank ligands and discover cryptic pockets.

**Article**: https://nostrumbiodiscovery.github.io/papers/Software/index.html#frag-pele-dynamic-ligand-growing-within-a-binding-site-a-novel-tool-for-hit-to-lead-drug-design

**Input** (further explained below):


    - protein-fragment_core pdb
    - Grown ligands sdf file

**Output**:

    - Ranking of the fragments
    - Ranked binding modes within each fragment

**Computational time**:

    - 30min/3h per fragment (depends on how many rotatable bonds the fragment has)

1. Complex Preparation
======================
   
Prepare the system with maestro (Protein Preparation Wizard)
and output a complex.pdb. The complex.pdb must contain protein-core_fragment in the desired initial configuration (Core fragment dock to the protein).

Make sure the ligand has:

 - Unique chain
 - No atomnames with spaces or single letter
 - Any residuename except UNK
 - Well defined aromatic bonds

2. Input Preparation
=====================
 
Prepare the input file ``input.yml``:

From sdf
+++++++++++++++++++++++++++++++++++++

..  code-block:: yaml

    frag_core: "/home/daniel/PR_core.pdb" #Protein-core complex
    frag_ligands: "/home/daniel/grown_ligands.sdf" #Full grown ligands
    resname: "LIG" #Ligand core residue name
    chain_core: "L" #Ligand core chain
    cpus: 48

For more optional flags please refer to `optative falgs <../../documentation/index.html>`_


3. Run simulation
====================


To run the system launch the simulation with the next command:

``python -m pele_platform.main input.yml``

4. Output
===============

Rank fragments:

``simulation_score.tsv``

Ranked top poses of each fragment:

``fragment_folder/top_results/``


