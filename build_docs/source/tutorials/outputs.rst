Understanding the output files
================================

Folder structure
-----------------

Each time you launch a **new simulation**, the Platform will generate a top level directory including all input and output
files.

Unless you provide a **custom folder name** using ``working_folder`` parameters, it will automatically enumerate one
based on the name of the perturbed ligand, e.g. ``LIG_Pele``, ``LIG_Pele_1``, etc. ensuring the old simulations will not
be overwritten.

.. code:: console

    .
    └── LIG_Pele
        ├── Data
        ├── DataLocal
        ├── Documents
        ├── input
        ├── output
        └── results
            ├── top_poses
            ├── plots
            └── clusters

Top level directory
-------------------

The top level directory contains the logs, configuration files and a copy of input.yaml created by the user.

    * ``LIG.log:`` log file
    * ``adaptive.conf``: control file for AdaptivePELE (`learn more <https://adaptivepele.github.io/AdaptivePELE/Examples.html#control-file-outline>`_)
    * ``input.yaml``: user's configuration file
    * ``pele.conf``: control file for PELE (`learn more <https://nostrumbiodiscovery.github.io/pele_docs/GeneralStructure/GeneralStructure.html>`_)
    * ``LIG.mae``: the rotamer library

Additionally, it contains several Data and Documents folder required by PELE++, which contain residue templates and
rotamer libraries. More details `here <https://nostrumbiodiscovery.github.io/pele_docs/molecularParameters.html>`_.

Input
-----

The input folded contains:

    * ``system_preprocessed.pdb``: original PDB file preprocessed by the Platform
    * ``input*.pdb``: input files with initial poses (depends on the package)
    * ``receptor.pdb``: PDB file with extracted protein
    * ``ligand.pdb``: PDB file with extracted ligand.

Output
------

Most importantly, this directory contains **raw PELE output**, i.e. report and trajectory files containing every step of the simulation, split
into several folders, depending on the number of iterations performed.

Additionally, you will find information about topologies and clustering performed by AdaptivePELE.

Results
-------

The results directory holds a more **user-friendly, curated output** comprising of three separate folders:

    * ``top_poses`` - top 100 lowest binding energy structures, unless the user specified a different number of top poses to be extracted
    * ``plots`` - plots of multiple metrics selected by the user (such as SASA or atom distance), providing an insight into the progress of the simulation
    * ``clusters`` - lowest binding energy cluster representatives (unless the user specified a different metric) and clustering plots.
