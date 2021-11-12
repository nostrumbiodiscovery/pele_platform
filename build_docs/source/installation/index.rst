Installation
===============


Conda
-----

.. code-block:: bash

    # Create a new conda environment from scratch:
    conda create -n pele_platform

    # Activate the new conda environment:
    conda activate pele_platform

    # Install Python 3.8:
    conda install -c conda-forge python=3.8

    # Install PELE Platform 1.6.1
    conda install -c nostrumbiodiscovery -c conda-forge pele_platform=1.6.1

    # Export environment variables
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
