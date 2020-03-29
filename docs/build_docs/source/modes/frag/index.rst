FragPELE Modes
######################



Summary of all FragPele growing methods


Standard
--------------

This method works by specifying the inicial receptor-core pdb and
an input file that contains information of the where and what to grow

- **frag_core**: Core of the molecule we want to add fragments to. Required parameter

- **frag_input**: Serie file with growing instructions. For more please refer here.

..  code-block:: yaml

    frag_core: "/home/daniel/PR_core.pdb"
    fag_input: "/home/daniel/serie_file.conf"






Sdf with full ligands
------------------------


This method works by specifying the inicial receptor-core pdb and
a sdf file with the full ligand. **All ligands must have molecule name**

- **frag_core**: Core of the molecule we want to add fragments to. Required parameter

- **frag_ligands**: Sdf with the aimed grown ligands

..  code-block:: yaml

    frag_core: "/home/daniel/PR_core.pdb"
    frag_ligands: "/home/daniel/grown_ligands.sdf"


AI
-----------

This method works by specifying the inicial receptor-core. Then a RNN 
generative model will grow fragments with as many atoms as iterations
the user set. i.e iterations: 3 will grow fragments up to 3 atoms all around
the molecule.

- **frag_core**: Core of the molecule we want to add fragments to. Required parameter

- **frag_ai**: Whether to use AI method or not

- **iterations**: Maximum number of atoms of your fragment

..  code-block:: yaml

    frag_core: "/home/daniel/PR_core.pdb"
    frag_ai: true 
    iterations: 7 
