Side chain perturbation parameters
----------------------------------

These are parameters that affect the side chain perturbation algorithm of
PELE.

List of PELE general parameters:

    1. `covalent_residue <#covalent-residue>`__
    2. `perturbation_trials <#perturbation-trials>`__
    3. `refinement_angle <#refinement-angle>`__

List of examples:

    - `Example 1 <#example-1>`__


covalent_residue
++++++++++++++++

    - Description: Sets chain and residue number of the residue to perturb
      with the side chain perturbation algorithm.
    - Type: ``Character``:``Integer``, chain and residue number
    - Default: ``None``

    .. note::
       This parameter activates the ``sidechainPerturbation`` package. Thus,
       it is a special case that do not need any further package to be chosen.

    .. seealso::
      `perturbation_trials <#perturbation-trials>`__,
      `refinement_angle <#refinement-angle>`__,
      `packages <../packages.html>`__,
      `Example 1 <#example-1>`__


perturbation_trials
+++++++++++++++++++

    - Description: Sets the number of trials for the side chain perturbation
      algorithm.
    - Type: ``Integer``
    - Default: ``20``

    .. seealso::
      `covalent_residue <#covalent-residue>`__,
      `refinement_angle <#refinement-angle>`__,
      `packages <../packages.html>`__,
      `Example 1 <#example-1>`__


refinement_angle
++++++++++++++++

    - Description: Sets the refinement angle which affects to the refinement
      stage. During the refinement, the side chain perturbation algorithm
      will perform a very soft perturbation, taking random angles between
      (-refinement_angle, +refinement_angle).
    - Type: ``Integer``
    - Default: ``10``

    .. seealso::
      `covalent_residue <#covalent-residue>`__,
      `perturbation_trials <#perturbation-trials>`__,
      `packages <../packages.html>`__,
      `Example 1 <#example-1>`__


Example 1
+++++++++

In this example we set a side chain perturbation simulation using 30
computation cores. This simulation will consist on perturbing the side
chains of the chosen residue. It can be a standard residue or non
standard like an amino acid covalently attached to a small molecule.
We also modify the default parameters of the side chain perturbation
algorithm, increasing the number of trials and the refinement angle
we can achieve more exploration.

..  code-block:: yaml

    # Required parameters
    system: 'system.pdb'
    chain: 'A'
    resname: 'LIG'

    # General parameters
    cpus: 30
    seed: 2021

    # Package selection
    covalent_residue: "A:273"

    # Side chain perturbation parameters
    perturbation_trials: 12
    refinement_angle: 15
