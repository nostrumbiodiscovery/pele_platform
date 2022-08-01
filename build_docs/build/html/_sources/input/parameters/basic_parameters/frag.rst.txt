fragPELE parameters
-------------------

These are parameters to set up a fragPELE simulation.

List of fragPELE parameters:

    1. `frag_core <#frag-core>`__
    2. `chain_core <#chain-core>`__
    3. `frag_ligands <#frag-ligands>`__
    4. `frag_library <#frag-library>`__
    5. `frag_core_atom <#frag-core-atom>`__
    6. `fragment_atom <#fragment-atom>`__
    7. `growing_steps <#growing-steps>`__
    8. `steps_in_gs <#steps-in-gs>`__
    9. `sampling_steps <#sampling-steps>`__
    10. `pele_control_file <#pele_control_file>`__

List of examples:

    - `Example 1 <#example-1>`__
    - `Example 2 <#example-2>`__
    - `Example 3 <#example-3>`__
    - `Example 4 <#example-4>`__
    - `Example 5 <#example-5>`__


frag_core
+++++++++

    - Description: Defines the path to PDB file containing the protein
      and docked scaffold.
    - Type: ``String``
    - Default: ``None``

    .. note::
       It is a mandatory parameter for fragPELE.

    .. note::
       The scaffold or core is the molecule on which the fragments will
       be inserted.

    .. seealso::
      `Example 1 <#example-1>`__,
      `Example 4 <#example-4>`__


chain_core
++++++++++

    - Description: Sets the unique chain id of the scaffold supplied
      in the ``frag_core`` structure.
    - Type: ``Character``
    - Default: ``None``

    .. note::
       It is a mandatory parameter for fragPELE.

    .. note::
       The scaffold, also referred to as core, is the molecule on which the
       fragments will be inserted.

    .. seealso::
      `frag_core <#frag-core>`__,
      `Example 1 <#example-1>`__,
      `Example 4 <#example-4>`__


frag_ligands
++++++++++++

    - Description: Defines the path to the SDF containing fully grown
      ligands. Fully grown means that each ligand must contain both
      scaffold and fragment already attached.
    - Type: ``String``
    - Default: ``None``

    .. warning::
       Each ligand present in the SDF must contain the scaffold defined
       in the ``frag_core`` and one fragment. If the substructure searcher
       is not able to find the scaffold in any of the fragments, the
       simulation will fail.

    .. note::
       There is an alternative method to run fragPELE. Instead of supplying
       an SDF containing fully grown ligands, we can supply a library of
       fragments with a parameter called ``frag_library``.

    .. note::
       Note that ``frag_ligands`` and ``frag_library`` parameters cannot be
       defined simultaneously.

    .. seealso::
      `frag_core <#frag-core>`__,
      `chain_core <#chain-core>`__,
      `frag_library <#frag-library>`__,
      `Example 1 <#example-1>`__,
      `Example 4 <#example-4>`__


frag_library
++++++++++++

    - Description: Defines the path to a folder containing fragment files.
      Fragments can be supplied as PDB or SDF but all of them must be
      placed into the same folder. Each fragment included into
      the ``frag_library`` directory will be inserted to the scaffold
      (or molecular core). The chemical bonding will take place between
      the atom of the scaffold selected with ``frag_core_atom`` and
      all non symmetric hydrogen atoms found in each fragment. We can
      also fix the hydrogen atom of the fragment we want to connect
      setting the parameter called ``fragment_atom``.
    - Type: ``String``
    - Default: ``None``

    .. note::
       There is an alternative method to run fragPELE. Instead of supplying
       a fragment library, we can supply an SDF containing fully grown
       ligands, with a parameter called ``frag_ligands``.

    .. note::
       Note that ``frag_ligands`` and ``frag_library`` parameters cannot be
       defined simultaneously.

    .. seealso::
      `frag_core <#frag-core>`__,
      `chain_core <#chain-core>`__,
      `frag_ligands <#frag-ligands>`__,
      `frag_core_atom <#frag-core-atom>`__,
      `fragment_atom <#fragment-atoms>`__,
      `Example 2 <#example-2>`__,
      `Example 3 <#example-3>`__


frag_core_atom
++++++++++++++

    - Description: Defines which is the atom of the scaffold the fragments
      must be connected to when using fragment libraries
      (``frag_library`` parameter).
    - Type: Two atoms, ``String``-``String``, where each string corresponds
      to the name of one atom of the scaffold.
        - First atom: heavy atom connected to the hydrogen atom that will
          be replaced with each fragment.
        - Second atom: hydrogen atom to replace.
    - Default: ``None``

    .. note::
       Note that this parameter only has an effect when a fragment library
       is supplied through the ``frag_library`` parameter.

    .. note::
       It is a mandatory parameter for fragPELE when a fragment library
       is supplied.

    .. seealso::
      `frag_library <#frag-library>`__,
      `Example 2 <#example-2>`__,
      `Example 3 <#example-3>`__


fragment_atom
+++++++++++++

    - Description: Defines which is the atom of each fragment the scaffold
      must be connected to when using fragment libraries
      (``frag_library`` parameter).
    - Type: One atom, ``String``, hydrogen atom to remove and replace
      with the scaffold.
    - Default: ``None``

    .. warning::
       When ``fragment_atom`` is specified, all fragments from the library
       must contain one hydrogen atom that matches with that name. Then,
       the connection to the scaffold will be applied through that position.
       This strategy requires a manual selection of each attachment atom
       and the assignment of the right PDB atom name to it.

    .. note::
       Note that this parameter only has an effect when a fragment library
       is supplied through the ``frag_library`` parameter.

    .. note::
       It is an optional parameter. When missing, bonding to the scaffold
       will take place through all asymmetric hydrogen atoms.

    .. seealso::
      `frag_library <#frag-library>`__,
      `Example 3 <#example-3>`__


growing_steps
+++++++++++++

    - Description: Sets the number of growing steps to apply during the
      growth of the fragment.
    - Type: ``Integer``
    - Default: ``6``

    .. note::
       Increasing the number of growing steps will smooth the alchemical
       change during the growth of the fragment but the simulation will
       become more expensive.

    .. seealso::
      `steps_in_gs <#steps-in-gs>`__,
      `sampling_steps <#sampling-steps>`__,
      `Example 3 <#example-3>`__


steps_in_gs
+++++++++++

    - Description: Sets the number of PELE steps to perform at each
      growing step.
    - Type: ``Integer``
    - Default: ``3``

    .. note::
       Increasing the number of growing steps will promote the conformational
       sampling and reallocation of the ligand and its neighboring side
       chains but the simulation will become more expensive.

    .. seealso::
      `growing_steps <#growing-steps>`__,
      `sampling_steps <#sampling-steps>`__,
      `Example 3 <#example-3>`__


sampling_steps
++++++++++++++

    - Description: Sets the number of PELE steps to perform during the
      final equilibration stage, which happens once the fragment is fully
      grown.
    - Type: ``Integer``
    - Default: ``20``

    .. note::
       Increasing the number of equilibration steps will promote the conformational
       sampling and reallocation of the ligand and its neighboring side
       chains but the simulation will become more expensive.

    .. seealso::
      `growing_steps <#growing-steps>`__,
      `steps_in_gs <#steps-in-gs>`__,
      `Example 3 <#example-3>`__


pele_control_file
+++++++++++++++++

    - Description: Sets a custom control file template for PELE that
      will replace the predefined template that fragPELE uses.
    - Type: ``str``
    - Default: ``None``

    .. note::
       The template must have certain parameters assigned through
       predetermined flags (marked with the dollar symbol: ``$``) so
       fragPELE can change them dynamically.
       Check an example of a template here:
       :download:`pele_template.conf <../../../../../_static/files/pele_template.conf>`

    .. seealso::
      `Example 5 <#example-5>`__


Example 1
+++++++++

In this example we set up a fragPELE simulation with 30 computation
cores. The goal is to take the initial structure supplied with the
``frag_core`` parameter and alchemically convert it to molecules
defined with the ``frag_ligands`` parameter.

..  code-block:: yaml

    # Required parameters
    frag_core: "complex_with_scaffold.pdb"
    chain_core: "L"
    resname: "LIG"

    # General parameters
    cpus: 30
    seed: 2021

    # fragPELE parameters
    frag_ligands: "fully_grown_ligands.sdf"


Example 2
+++++++++

In this example we set up a fragPELE simulation with 30 computation
cores. The goal is to take the initial structure supplied with the
``frag_core`` parameter and alchemically attach all fragments defined
in the library files from the path set by the ``frag_library`` parameter.
We must also specify the atom of the scaffold where fragments need to
be inserted using ``frag_core_atom`` parameter. In this case, we
attach fragments through a hydrogen atom called **H6** that is connected
to a carbon atom with name **C6**. Fragments will be connected to
this position through all asymmetric hydrogen atoms.

..  code-block:: yaml

    # Required parameters
    frag_core: "complex_with_scaffold.pdb"
    chain_core: "L"
    resname: "LIG"

    # General parameters
    cpus: 30
    seed: 2021

    # fragPELE parameters
    frag_library: "path/to/frag/libraries"
    frag_core_atom: "C6-H6"


Example 3
+++++++++

In this example we set up a fragPELE simulation with 30 computation
cores. The goal is to take the initial structure supplied with the
``frag_core`` parameter and alchemically attach all fragments defined
in the library files from the path set by the ``frag_library`` parameter.
We must also specify the atom of the scaffold where fragments need to
be inserted using ``frag_core_atom`` parameter. In this case, we
attach fragments through a hydrogen atom called **H6** that is connected
to a carbon atom with name **C6**. Since we also supply the ``fragment_atom``
parameter, fragments will be connected to atom **C6** from scaffold
through the hydrogen atom called **HGRW**.

..  code-block:: yaml

    # Required parameters
    frag_core: "complex_with_scaffold.pdb"
    chain_core: "L"
    resname: "LIG"

    # General parameters
    cpus: 30
    seed: 2021

    # fragPELE parameters
    frag_library: "path/to/frag/libraries"
    frag_core_atom: "C6-H6"
    fragment_atom: "HGRW"


Example 4
+++++++++

In this example we set up a fragPELE simulation with 30 computation
cores. The goal is to take the initial structure supplied with the
``frag_core`` parameter and alchemically convert it to molecules
defined with the ``frag_ligands`` parameter. Besides, we are
significantly increasing the length of the alchemical growth because
we ask for more growing steps (``growing_steps``) and more PELE steps
per growing step (``steps_in_gs``). On the other hand, we reduce the
length of the final equilibration (``sampling_steps``).

..  code-block:: yaml

    # Required parameters
    frag_core: "complex_with_scaffold.pdb"
    chain_core: "L"
    resname: "LIG"

    # General parameters
    cpus: 30
    seed: 2021

    # fragPELE parameters
    frag_ligands: "fully_grown_ligands.sdf"
    growing_steps: 10
    steps_in_gs: 5
    sampling_steps: 10


Example 5
+++++++++

In this example we set up a fragPELE simulation with 48 computation
cores. The goal is to take the initial structure supplied with the
``frag_core`` parameter and alchemically convert it to molecules
defined with the ``frag_ligands`` parameter. Besides, we ask
to use peleffy along with the Open Force Field parameters for
hetero molecules with ``use_peleffy`` and ``forcefield``
parameters. Finally, we replace fragPELE's default control file
template with another template that we sett with
``pele_control_file``.

..  code-block:: yaml

    # Required parameters
    frag_core: "complex_with_scaffold.pdb"
    chain_core: "L"
    resname: "LIG"

    # General parameters
    cpus: 48
    seed: 2022

    # fragPELE parameters
    frag_ligands: "fully_grown_ligands.sdf"
    use_peleffy: True
    forcefield: 'openff-2.0.0'
    pele_control_file: "pele_template.conf"
