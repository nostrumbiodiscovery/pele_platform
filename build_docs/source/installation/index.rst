Installation
===============


Conda
-----

.. code-block:: bash

    conda install -c nostrumbiodiscovery -c conda-forge pele_platform
    
    export PELE=/path/to/PELE-1.X/

    export SCHRODINGER=/path/to/schoringer/20XX/

    export PELE_LICENSE=/path/to/pele/licenses/

For more information on setting up your environment and licenses, see `Environment Variables <../environment/index.html>`_


Singularity
-----------

We have several ways available to run pele_platform from Singularity containers. This method is ideal for running pele_platform and PELE software on HPC systems.

If you are interested in learning more about it, please contact it@nostrumbiodiscovery.com.


Test installation
--------------------

.. code-block:: bash

    git clone https://github.com/NostrumBioDiscovery/pele_platform.git

    cd pele_platform/tests

    python -m pele_platform.main ../pele_platform/Examples/induced_fit/input_fast.yaml
