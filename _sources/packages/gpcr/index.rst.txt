Prepare your own GPCR orthosteric site simulation
#########################################################################

This simulation aims to find the binding path
of a small molecule to the orthosteric site of a GPCR.

**Input** (further explained below):

    - protein-ligand pdb
    - Orthosteric site atom
    - Initial site atom

**Output** (further explained below):

    - ranked binding modes of the chosen small molecule
      with the desired protein

**Computational Time**: 3h

1. Complex Preparation
======================
   
Prepare the system with maestro (Protein Preparation Wizard)
and output a complex.pdb. The complex.pdb must contain the protein-ligand.  The ligand can be place anywhere as it will be automatically placed all around the specified initial site atom (more below).

2. Input Preparation
=====================

This simulation will initially place the ligand all over the initial_site. Then a slight
bias will be performed against SASA and a small box containing the initial and orthosteric
site will be set for simulation. High constraints will be placed on the CA to avoid 
the initial structure to collapse because of the lack of the membrane.

Prepare the input file ``input.yml``:

..  code-block:: yaml

    system: "complex.pdb" #Ligand-protein complex
    resname: "LIG" #Ligand residue name in system
    chain: "L" #Ligand chain name in syste,
    seed: 1234
    gpcr_orth: true #Set defaults for gpcr simulation
    orthosteric_site: "A:114:CG" #Atom in the orthosteric site
    initial_site: "A:310:CD" #Atom in the initial site
    cpus: 50

For more optional flags please refer to `optative falgs <../../documentation/index.html>`_


3. Run simulation
====================

To run the system launch the simulation with the next command:

``python -m pele_platform.main input.yml``

4. Output
=================

Best ranked clusters:

``working_folder/results/clusters``

Best ranked poses:

``working_folder/results/BestStructs/``
