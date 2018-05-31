===============
Getting Started
===============

.. toctree::
   :maxdepth: 2

These MSM Pele cookbook shows how to use basics of the platform and walk you
through your first steps.

Launch your first MSM job:
----------------------------

**Previous Requisites**

- ligand-receptor complex or both separate.

**Launch MSM_Pele with OPLS charges**::

    python2.X main.py PDB ligand_resname ligand_chain --cpus X


**Launch MSM_Pele with QM charges**::

    python2.X main.py receptor_PDB ligand_resname ligand_chain --mae_lig ligand.mae --cpus X

