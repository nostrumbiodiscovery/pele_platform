Environment variables
======================

Mandatory variables
-------------------

They must be defined for the pele_platform to work correctly:
	- PELE: path to PELE installation
	- SCHRODINGER: path to Schrodinger installation

	Examples:

	.. code-block:: bash

		export PELE=/path/to/PELE-1.X/
		export SCHRODINGER=/path/to/schoringer/20XX/

Optional variables
------------------

PELE_EXEC:
++++++++++

    The pele executable can be configured in several ways (from highest to lowest precedence):

      1. Can be specified by input.yaml flag:

      .. code-block:: yaml

        pele_exec: "/home/pele/bin/Pele_mpi"

      2. Can be specified by environment variable:

      .. code-block:: bash

        export PELE_EXEC=/home/pele/bin/Pele_mpi

      3. If no path is specified, the pele executable must be /path/to/PELE-1.X/bin/Pele_mpi (according to the path specified by the environment variable PELE)

PELE_DATA:
++++++++++

    The pele data folder can be configured in several ways (from highest to lowest precedence):

      1. Can be specified by input.yaml flag:

      .. code-block:: yaml

        pele_data: "/home/pele/Data/"

      2. Can be specified by environment variable:

      .. code-block:: bash

        export PELE_DATA=/home/pele/Data/

      3. If no path is specified, the pele data folder must be /path/to/PELE-1.X/Data (according to the path specified by the environment variable PELE)

PELE_DOCUMENTS:
+++++++++++++++

    The pele documents folder can be configured in several ways (from highest to lowest precedence):

      1. Can be specified by input.yaml flag:

      .. code-block:: yaml

        pele_documents: "/home/pele/Documents/"

      2. Can be specified by environment variable:

      .. code-block:: bash

        export PELE_DOCUMENTS=/home/pele/Documents/

      3. If no path is specified, the pele documents folder must be /path/to/PELE-1.X/Documents (according to the path specified by the environment variable PELE)

PELE_LICENSE:
+++++++++++++

	Licenses can be configured in several ways (from highest to lowest precedence):

	  1. Can be specified by input.yaml flag:

	  .. code-block:: yaml

	  	license: /path/to/licenses/folder/

	  2. Can be specified by environment variable:

	  .. code-block:: bash

	  	export PELE_LICENSE=/path/to/licenses/folder/

      3. If no path is specified, licenses must be under /path/to/PELE-1.X/licenses/ (according to the path specified by the environment variable PELE)


SINGULARITY_EXEC:
+++++++++++++++++

	Path of the singularity container that contains the PELE executable. It can be configured in two ways:

	  1. Can be specified by input.yaml flag:

	  .. code-block:: yaml

	  	singularity_exec: /path/to/singularity/container/


	  2. Can be specified by environment variable:

	  .. code-block:: bash

		export SINGULARITY_EXEC=/path/to/singularity/container/
