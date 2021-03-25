Installation
===============

Conda (recommended)
-----------------------

.. code-block:: bash

    conda install -c nostrumbiodiscovery -c conda-forge -c anaconda pele_platform==1.5.1
    
    export PELE=/path/to/PELE-1.X/

    export SCHRODINGER=/path/to/schoringer/20XX/

    Licenses must be under /path/to/PELE-1.X/licenses/. Otherwise, they can be
    specifed via input.yaml flag. i.e license: /path/to/licenses/folder/


Singularity
----------------

.. code-block:: bash

   git clone https://github.com/NostrumBioDiscovery/pele_platform.git

   cd pele_platform/singularity/
   
   bash create_pele_image.sh (You will need acces to PELE dockerhub)
   
   cd test/
   
   # Change the licenses and folders inside singularity.sl (IT section)
   
   bash singularity.sl  # If not queue system
   
   sbatch singularity.sl  # If slurm queue system



Pypi
------

.. code-block:: bash

    pip install numpy cython

    pip install pele_platform==1.5.1

    export PELE=/path/to/PELE-1.X/

    export SCHRODINGER=/path/to/schrodinger/20XX/
    
    conda install rdkit (or compile from source code)

    Licenses must be under /path/to/PELE-1.X/licenses/. Otherwise, they can be
    specified via input.yaml flag, i.e license: /path/to/licenses/folder/


Last stable release from source code
--------------------------------------------

.. code-block:: bash

    git clone https://github.com/NostrumBioDiscovery/pele_platform.git
    
    cd pele_platform
    
    pip install .
    
    export PELE=/path/to/PELE-1.X/

    export SCHRODINGER=/path/to/schoringer/20XX/

    conda install rdkit (if you want to have the possibility to build constraints by SMILES)

    Licenses must be under /path/to/PELE-1.X/licenses/. Otherwise, they can be
    specified via input.yaml flag, i.e license: /path/to/licenses/folder/



Test installation
----------------

.. code-block:: bash

    git clone https://github.com/NostrumBioDiscovery/pele_platform.git

    cd pele_platform/tests

    python -m pele_platform.main ../pele_platform/Examples/induced_fit/input_fast.yaml
