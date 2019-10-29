Installation
###############

Conda (recomended)
--------------------

on going


Pypi
------

.. code-block:: bash

    pip install numpy cython

    pip install pele_platform


Source Code
-------------

.. code-block:: bash

    git clone https://github.com/NostrumBioDiscovery/pele_platform.git
    
    cd pele_platform
    
    vim pele_platform/constants/constants.py #(change paths under else statement)
    
    pip install pele_platform



Test it works
----------------

.. code-block:: bash

    cd pele_platform/tests

    pytest
