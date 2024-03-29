================================================
Large protein conformational moves with PCA-PELE
================================================

Introduction
-------------------

This package aims to analyse and enhance the receptor movement by performing the principal component analysis (PCA)
on three or more trajectory snapshots, which is especially useful when the protein undergoes big conformational changes.
The best binding modes of the simulated ligand are ranked and extracted for user inspection.

Inputs
+++++++++

    - protein-ligand PDB file
    - 3 identical PDB files (only changing coordinates) which define a motion
    - YAML file with parameters

Default parameters
+++++++++++++++++++

Exact parameters depend on the selected simulation type, the PCA only enhances the receptor movements.

**Computational time**: 5h 

1. Complex Preparation
---------------------------
   
Prepare the system with Maestro (Protein Preparation Wizard) and output a ``complex.pdb``. The file must contain the
protein-ligand in the desired initial configuration.

Make sure the ligand has:

 - unique chain ID
 - no atom names with spaces or single letter
 - any residue name except ``UNK``

2. Input Preparation
---------------------------

Prepare the input file ``input.yml``:

..  code-block:: yaml

    system: 'complex.pdb' # protein-ligand complex
    chain: 'L' # chain name of the ligand
    resname: 'LIG' # residue name of the ligand
    seed: 12345

    # Calculate PCA by giving 3 pdb snapshots
    # defining a motion
    pca_traj:
    - "snap1.pdb"
    - "snap2.pdb"
    - "snap3.pdb"
    # Do not constraint structure
    remove_constraints: true
    # Big sampling of the receptor while making binding
    spawning: independent
    out_in: true #You can change the method to any other (induced_fit, rescoring...)

For more optional flags please refer to `optional flags <../../input/yaml.html>`_.

3. Run simulation
--------------------

To run the system launch the simulation with the following command:

``python -m pele_platform.main input.yml``

4. Output
-----------

Clusters
++++++++++

Upon completion of the simulation, all trajectories are clustered based on ligand heavy atom coordinates. Then, a cluster representative with the best binding energy (or metric of your choice) is selected.
Ranked cluster representatives can be found in:

``working_folder/results/clusters``

Best snapshots
+++++++++++++++++

In addition, top 100 structures with the best binding energy (or metric of your choice) are retrieved. This is done to ensure the clustering algorithm did not skip any valuable results. They are stored in:

``working_folder/results/top_poses``
