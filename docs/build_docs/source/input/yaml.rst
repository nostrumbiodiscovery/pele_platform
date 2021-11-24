=====================
Input yaml parameters
=====================

This section contains the full list of parameters that can be inserted to ``input.yaml``
in order to change the behavior of the different algorithms of PELE.


Required parameters
-------------------

These parameters must always be present in the ``input.yaml``. They define the basic
settings that the Platform needs in order to run.

Any PELE simulation needs, at least:

    - ``system``: the path to the PDB file that contains the initial structure
      for the simulation. It can also take several PDB files with the asterisk
      wildcard.
    - ``chain``: the chain in the input PDB file that corresponds to the small
      molecule to perturb.
    - ``resname``: the residue name in the input PDB file that corresponds to
      the small molecule to perturb.

To use a single input PDB file:

..  code-block:: yaml

  system: 'system.pdb'
  chain: 'L'
  resname: 'LIG'


To take all PDB files inside a folder:

..  code-block:: yaml

  system: 'PDBs/*.pdb'
  chain: 'L'
  resname: 'LIG'


.. warning::
  Note that fragPELE uses different parameters than the rest of packages
  of PELE. Please, refer to fragPELE page to see the parameters that are required
  in that specific case.


Packages
--------

To define the type of simulation to run, a package needs to be selected in
the ``input.yaml``. There are several packages available in the Platform.
Some of them can be activated by setting specific parameters to
true. They are listed below:

    - Local explorations:

        - ``induced_fit_fast``: fast induced-fit docking protocol.
          More information `here <../packages/induced_fit/index.html>`_.
        - ``induced_fit_long``: more expensive induced-fit docking protocol.
          More information `here <../packages/induced_fit/index.html>`_.
        - ``rescoring``: local conformation refinement.
          More information `here <../packages/rescoring/index.html>`_.

    - Global explorations:
        - ``out_in``: ligand migration from the bulk of solvent to a protein cavity.
          More information `here <../packages/migration/binding.html>`_.
        - ``in_out``: ligand migration from a protein cavity to the bulk of solvent.
          More information `here <../packages/migration/unbinding.html>`_.
        - ``site_finder``: global search of ligand binding sites.
          More information `here <../packages/site_finder/index.html>`_.

Other packages can be activated when specific parameters are set, like:

    - fragPELE: our fragment growing protocol has its own parameters.
      More information `here <../packages/frag/index.html>`_.
    - sidechainPerturbation: when ``covalent_residue`` is set, the side chain
      perturbation algorithm will run to perturb covalently bound non standard
      residues. More information `here <../packages/site_finder/index.html>`_.


Basic parameters
----------------

Find below a list of all

.. toctree::
   :maxdepth: 1

   basic_parameters/general.rst
   basic_parameters/pele.rst
   basic_parameters/adaptive.rst



fragPELE parameters
-------------------

.. toctree::
   :maxdepth: 2
   :caption: fragPELE parameters

   frag/index.rst


Advanced parameters
-------------------

.. toctree::
   :maxdepth: 2
   all_packages/index.rst
   :caption: Advanced parameters
