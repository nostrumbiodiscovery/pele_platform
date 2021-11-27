"""""""""""""""""
Package selection
"""""""""""""""""

To define the type of simulation to run, a package needs to be selected in
the ``input.yaml``. There are several packages available in the Platform.
Some of them can be activated by setting specific parameters to
True. They are listed below:

    - Local explorations:

        - ``induced_fit_fast``: fast induced-fit docking protocol.
          More information `here <../../packages/induced_fit/index.html>`__.
        - ``induced_fit_long``: more expensive induced-fit docking protocol.
          More information `here <../../packages/induced_fit/index.html>`__.
        - ``rescoring``: local conformation refinement.
          More information `here <../../packages/rescoring/index.html>`__.

    - Global explorations:
        - ``out_in``: ligand migration from the bulk of solvent to a protein cavity.
          More information `here <../../packages/migration/binding.html>`__.
        - ``in_out``: ligand migration from a protein cavity to the bulk of solvent.
          More information `here <../../packages/migration/unbinding.html>`__.
        - ``site_finder``: global search of ligand binding sites.
          More information `here <../../packages/site_finder/index.html>`__.

Other packages can be activated when specific parameters are set, like:

    - fragPELE: our fragment growing protocol has its own parameters.
      More information `here <basic_parameters/frag.html>`__.
    - sidechainPerturbation: when ``covalent_residue`` is set, the side chain
      perturbation algorithm will run to perturb covalently bound non standard
      residues.
      More information `here <basic_parameters/sidechain_perturbation.html>`__.
