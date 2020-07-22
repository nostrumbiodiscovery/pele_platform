FragPELE Modes
######################



Summary of all FragPele growing methods



HT Mode
------------------------


This method works by specifying the inicial receptor-core pdb and
a sdf file with the full ligand. **All ligands must have molecule name**


- **frag_core**: Core of the molecule we want to add fragments to. Required parameter

- **frag_ligands**: Sdf with the aimed grown ligands

- **resname**: Residue name of the frag_core ligand
 
- **cpus**: Cpus to use. Default: 48

- **chain_core**: Chain of the ligand core. Default: L

..  code-block:: yaml

    frag_core: "/home/daniel/PR_core.pdb"
    frag_ligands: "/home/daniel/grown_ligands.sdf"
    cpus: 48
    chain_core: "L"




LT Mode
--------------

This method works by specifying the inicial receptor-core pdb and
an input file that contains information of the where and what fragment to grow.
**The fragment need to be the chain L at the moment**

- **frag_core**: Core of the molecule we want to add fragments to. Required parameter

- **frag_input**: Serie file with growing instructions. For more please refer here.

- **cpus**: Cpus to use. Default: 48

- **chain_core**: Chain of the ligand core. Default: L

..  code-block:: yaml

    frag_core: "/home/daniel/PR_core.pdb"
    fag_input: "/home/daniel/serie_file.conf"
    cpus: 48
    chain_core: "L"
