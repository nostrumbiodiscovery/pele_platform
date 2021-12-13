"""""""""""""""""""
Required parameters
"""""""""""""""""""

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
  in this specific case.
