Prepare your own FragPELE simulation
######################################

FrAG, a new tool for in silico hit-to-lead drug design, capable of performing HT-fragment growing into a core while exploring the protein-ligand conformational space

1. Complex Preparation
======================
   
Prepare the system with maestro (Protein Preparation Wizard)
and output a complex.pdb. The complex.pdb must contain protein-ligand_to_be_grown in the desired initial configuration.

Make sure the ligand has:

 - Unique chain
 - No atomnames with spaces or single letter
 - Any residuename except UNK
 - Well defined aromatic bonds

2. Input Preparation
=====================
 
Prepare the input file ``input.yml``:

To run different modes prepare different control files.
For more explanation on the modes please refer to `here <../../modes/frag/index.html>`__


From sdf
+++++++++++++++++++++++++++++++++++++

..  code-block:: yaml

    frag_core: "/home/daniel/PR_core.pdb"
    frag_ligands: "/home/daniel/grown_ligands.sdf"
    resname: "LIG"
    cpus: 48
    chain_core: "L"

From serie file
+++++++++++++++++++++

..  code-block:: yaml

    frag_core: "../pele_platform/Examples/frag/1w7h_preparation_structure_2w.pdb"
    frag_input: "../pele_platform/Examples/frag/sequential.conf"
    resname: "LIG"
    cpus: 48
    chain_core: "L"

For more about the serie file please refer to: https://carlesperez94.github.io/frag_pele/first_steps.html


3. Run simulation
====================


To run the system launch the simulation with the next command:

``python -m pele_platform.main input.yml``

