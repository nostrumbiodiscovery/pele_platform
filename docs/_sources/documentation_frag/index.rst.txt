FragPELE Parameters
######################

Compulsory flags
--------------------

Frag PELE grows an atom onto a core in N growing steps while moving protein and ligand.
Afterwards a final sampling simulation is run to fully explore the ligand-protein conformational space.

- **frag_core**: Core of the molecule we want to add fragments to. Required parameter

- **fag_input**: Frag pele input. For more please check here

..  code-block:: yaml

    frag_core: "/home/daniel/PR_core.pdb"
    fag_input: "/home/daniel/serie_file.conf"

Optative flags
-------------------

- **growing_steps**: Number of steps to grow the fragment with.

- **steps_in_gs**: Number of pele steps within each growing step

- **sampling_steps**: Number of pele steps in the final sampling simulation

- **protocol**: Type of protocol. options = [HT, ES]. For more info please refere here.

- **cpus**: Cpus to use


..  code-block:: yaml

    growing_steps: 6
    steps_in_gs: 6
    sampling_steps: 20
    protocol: HT
    cpus: 24

