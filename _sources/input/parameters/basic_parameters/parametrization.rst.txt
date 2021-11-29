Residue parametrization
-----------------------

These are parameters to manage the parametrization of residues. For a complete
guide on the parameters that PELE needs in order to run successfully, check
the `Custom template generation <../../templates.html>`__ documentation.

List of ligand parameters:

    1. `use_peleffy <#use-peleffy>`__
    2. `templates <#templates>`__
    3. `rotamers <#rotamers>`__
    4. `solvent_template <#solvent-template>`__

List of examples:

    - `Example 1 <#example-1>`__
    - `Example 2 <#example-2>`__
    - `Example 3 <#example-3>`__
    - `Example 4 <#example-4>`__


use_peleffy
+++++++++++

    - Description: Activates peleffy instead of the default parameters builder.
    - Type: ``Boolean``
    - Default: ``False``

    .. warning::
       Although it is still not the default option, peleffy will soon be
       the default option to parametrize non standard residues.

    .. note::
       Note that the usage of any OpenFF force field requires activating
       peleffy.

    .. note::
       To get more information about peleffy you can check
       its `specific documentation  <https://martimunicoy.github.io/peleffy/index.html>`__
        or one of our `tutorials <../../../tutorials/peleffy.html>`__.

    .. seealso::
      `peleffy tutorial <../../../tutorials/peleffy.html>`__,
      `forcefield <pele.html#forcefield>`__,
      `Example 1 <#example-1>`__,
      `Example 2 <#example-2>`__,
      `Example 4 <#example-4>`__


templates
+++++++++

    - Description: Sets a list of external templates. Templates are
      files containing the parameters of the force field corresponding
      to a specific residue. Check `this page <../../templates.html>`__
      to get further information.
    - Type: list of ``String``
    - Default: ``None``

    .. warning::
       Double check that all atoms in the template match with those in the
       input structure of PELE for all involved residues. In case that
       PELE finds any mismatch between these files, it will raise an
       exception with information about any non matching atom that is found.
       To solve these problems, read carefully PELE messages and verify
       the parametrization process.

    .. note::
       Any residue parametrized in one of the external templates supplied
       with this flag will not be parametrized with our internal
       parametrization workflow.

    .. seealso::
      `rotamers <#rotamers>`__,
      `solvent_template <#solvent-template>`__,
      `Example 3 <#example-3>`__,
      `Example 4 <#example-4>`__


rotamers
++++++++

    - Description: Sets a list of rotamer libraries. Rotamer libraries are
      files containing the list of rotatable bonds of a residue arranged
      by rotatable branches. Check `this page <../../templates.html>`__
      to get further information.
    - Type: list of ``String``
    - Default: ``None``

    .. warning::
       When the template of a residue is externally supplied
       with ``templates`` option, the internal parametrization workflow
       skips the preparation of that residue. In order to account for
       its rotamers, you must also supply them through an external file as
       they will not be automatically generated.

    .. seealso::
      `templates <#templates>`__,
      `solvent_template <#solvent-template>`__,
      `Example 3 <#example-3>`__,
      `Example 4 <#example-4>`__


solvent_template
++++++++++++++++

    - Description: Sets a list of solvent templates. Solvent templates are
      files containing the solvent parameters of a residue.
      Check `this page <../../templates.html>`__
      to get further information.
    - Type: list of ``String``
    - Default: ``None``

    .. warning::
       When the template of a residue is externally supplied
       with ``templates`` option, the internal parametrization workflow
       skips the preparation of that residue. In order to correctly set
       its solvent parameters, you must also supply them through an external
       file as they will not be automatically generated.

    .. note::
       Solvent templates are only required when the OBC solvent model is used.

    .. seealso::
      `templates <#templates>`__,
      `rotamers <#rotamers>`__,
      `solvent <pele.html#templates>`__,
      `Example 3 <#example-3>`__,
      `Example 4 <#example-4>`__


Example 1
+++++++++

In this example we set an induced fit docking simulation with 30 computation
cores. We also activate peleffy and use the latest OpenFF force field
to parametrize non standard residues.

..  code-block:: yaml

    # Required parameters
    system: 'system.pdb'
    chain: 'L'
    resname: 'LIG'

    # General parameters
    cpus: 30
    seed: 2021

    # Package selection
    induced_fit_fast: True

    # PELE parameters
    forcefield: "openff-2.0.0"

    # Parametrization parameters
    use_peleffy: True


Example 2
+++++++++

In this example we set an induced fit docking simulation with 30 computation
cores. We also activate peleffy and use the OPLS2005 force field
to parametrize non standard residues.

..  code-block:: yaml

    # Required parameters
    system: 'system.pdb'
    chain: 'L'
    resname: 'LIG'

    # General parameters
    cpus: 30
    seed: 2021

    # Package selection
    induced_fit_fast: True

    # PELE parameters
    forcefield: "OPLS2005"

    # Parametrization parameters
    use_peleffy: True


Example 3
+++++++++

In this example we set an induced fit docking simulation with 30 computation
cores. We are using the default force field, which is OPLS2005. However,
in this case, we supply custom templates for our ligand. These files can
be generated externally with peleffy.

..  code-block:: yaml

    # Required parameters
    system: 'system.pdb'
    chain: 'L'
    resname: 'LIG'

    # General parameters
    cpus: 30
    seed: 2021

    # Package selection
    induced_fit_fast: True

    # Parametrization parameters
    templates:
        - "simulation_files/ligz"
    rotamers:
        - "simulation_files/LIG.rot.assign"


Example 4
+++++++++

In this example we set an induced fit docking simulation with 30 computation
cores. We are using the latest OpenFF force field, which works with the OBC
solvent model. So, in this case, if we supply custom templates for our ligand
we also need to supply the solvent template. These files can be generated
externally with peleffy.

..  code-block:: yaml

    # Required parameters
    system: 'system.pdb'
    chain: 'L'
    resname: 'LIG'

    # General parameters
    cpus: 30
    seed: 2021

    # Package selection
    induced_fit_fast: True

    # PELE parameters
    forcefield: "openff-2.0.0"

    # Parametrization parameters
    use_peleffy: True
    templates:
      - "simulation_files/ligz"
    rotamers:
      - "simulation_files/LIG.rot.assign"
    solvent_template:
      - "simulation_files/ligandParams.txt"
