Installation
###############

Conda (recomended)
--------------------

.. code-block:: bash

    conda install numpy cython

    conda install -c nostrumbiodiscovery -c conda-forge -c anaconda pele_platform

    vim /path/to/installation/pele_platform/constants/constants.py #(change paths under else statement)


Pypi
------

.. code-block:: bash

    pip install numpy cython

    pip install pele_platform
    
    vim /path/to/installation/pele_platform/constants/constants.py #(change paths under else statement)


Last stable release from source code
--------------------------------------------

.. code-block:: bash

    git clone https://github.com/NostrumBioDiscovery/pele_platform.git
    
    cd pele_platform
    
    pip install .
    
    vim /path/to/installation/pele_platform/constants/constants.py #(change paths under else statement)


Latest devel version
----------------------

.. code-block:: bash

    git clone https://github.com/NostrumBioDiscovery/pele_platform.git
    
    cd pele_platform

    git checkout devel
    
    pip install .

    vim /path/to/installation/pele_platform/constants/constants.py #(change paths under else statement)
    
Test it works
----------------

.. code-block:: bash

    cd pele_platform/tests

    pytest
