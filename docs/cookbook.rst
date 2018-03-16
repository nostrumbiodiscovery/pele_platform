===============
Getting Started
===============

.. toctree::
   :maxdepth: 2

Launch MSM_Pele with OPLS charges
-----------------------------------

The prerequisites to run MSM_Pele is:

    - PDB of the complex ligand+receptor

The launch commands is:
    
    - python2.X main.py PDB --cpus X


Launch MSM_Pele with QM charges
--------------------------------

The prerequisites to run MSM_Pele with keeping the QM charges of the ligand are:

    - PDB of the receptor
    - Maestro file of the ligand with the QM charges

The launch command is:

    - python2.X main.py receptor_PDB --mae_lig ligand.mae --cpus X

