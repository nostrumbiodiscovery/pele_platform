Installation
###############

Conda (recomended)
--------------------

.. code-block:: bash

    conda install numpy cython

    conda install -c nostrumbiodiscovery -c conda-forge -c anaconda pele_platform
    
    export PELE=/path/to/PELE-1.X/

    export SCHRODINGER=/path/to/schoringer/20XX/



Pypi
------

.. code-block:: bash

    pip install numpy cython

    pip install pele_platform

    export PELE=/path/to/PELE-1.X/

    export SCHRODINGER=/path/to/schoringer/20XX/
    
    conda install rdkit (if you want to have the possibility to build constraints by SMILES)


Last stable release from source code
--------------------------------------------

.. code-block:: bash

    git clone https://github.com/NostrumBioDiscovery/pele_platform.git
    
    cd pele_platform
    
    pip install .
    
    export PELE=/path/to/PELE-1.X/

    export SCHRODINGER=/path/to/schoringer/20XX/

    conda install rdkit (if you want to have the possibility to build constraints by SMILES)


Latest devel version
----------------------

.. code-block:: bash

    git clone https://github.com/NostrumBioDiscovery/pele_platform.git
    
    cd pele_platform

    git checkout devel
    
    pip install .

    export PELE=/path/to/PELE-1.X/

    export SCHRODINGER=/path/to/schoringer/20XX/


Test it works
----------------

.. code-block:: bash

    cd pele_platform/tests

    pytest
