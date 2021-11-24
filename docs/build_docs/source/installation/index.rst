============
Installation
============


PELE Platform is shipped as a Conda package that can be installed using ``conda-forge``
and ``nostrumbiodiscovery`` channels from Anaconda.


Requirements
------------

    - Miniconda3 or Conda package
    - Schrödinger Suite (2017-1 or later)
    - PELE++ code (optional, only required to run PELE but not to run PELE
      Plotter or the Analysis package)


OS support
----------
Currently, the platform only supports Linux and MacOS operating systems.


Conda installation
------------------

The most straightforward way to install the Platform is with Conda. The preferable
installation strategy is to generate a new Conda environment to hold the source code
and its dependencies. In order to generate a Conda environment and install the
Platform, please, follow the steps below:

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

    git clone https://github.com/nostrumbiodiscovery/pele_platform.git

    cd pele_platform/tests

    conda activate pele_platform

    python -m pele_platform.main ../pele_platform/Examples/induced_fit/input_fast.yaml
