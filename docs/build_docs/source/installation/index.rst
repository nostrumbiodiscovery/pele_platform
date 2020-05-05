Installation
###############

Conda (recommended)
-----------------------

.. code-block:: bash

    conda install -c nostrumbiodiscovery -c conda-forge -c anaconda pele_platform==1.4.4
    
    export PELE=/path/to/PELE-1.X/

    export SCHRODINGER=/path/to/schoringer/20XX/

    Licenses must be under /path/to/PELE-1.X/licenses/. Other non-default paths can be
    specifed via input.yaml flag. i.e license: /path/to/licenses/folder/



Pypi
------

.. code-block:: bash

    pip install numpy cython

    pip install pele_platform==1.4.4

    export PELE=/path/to/PELE-1.X/

    export SCHRODINGER=/path/to/schoringer/20XX/
    
    conda install rdkit (Or compile from source code)

    Licenses must be under /path/to/PELE-1.X/licenses/. Other non-default paths can be
    specifed via input.yaml flag. i.e license: /path/to/licenses/folder/


Last stable release from source code
--------------------------------------------

.. code-block:: bash

    git clone https://github.com/NostrumBioDiscovery/pele_platform.git
    
    cd pele_platform
    
    pip install .
    
    export PELE=/path/to/PELE-1.X/

    export SCHRODINGER=/path/to/schoringer/20XX/

    conda install rdkit (if you want to have the possibility to build constraints by SMILES)

    Licenses must be under /path/to/PELE-1.X/licenses/. Other non-default paths can be
    specifed via input.yaml flag. i.e license: /path/to/licenses/folder/



Test it works
----------------

.. code-block:: bash

    cd pele_platform/tests

    pytest test_*
