Environment variables
======================

Mandatory variables, they must be defined for the pele_platform works correctly:
	- PELE: path to PELE installation
	- SCHRODINGER: path to Schrodinger installation

	Examples:

	.. code-block:: bash

		export PELE=/path/to/PELE-1.X/
		export SCHRODINGER=/path/to/schoringer/20XX/

Optional variables:
    - PELE_LICENSE:

    	Licenses can be configured in several ways:

    	  1. Can be specified by input.yaml flag: 

    	  .. code-block:: yaml

    	  	license: /path/to/licenses/folder/

    	  2. Can be specified by environment variable:

    	  .. code-block:: bash

    	  	export PELE_LICENSE=/path/to/licenses/folder/

          3. If no path is specified, licenses must be under /path/to/PELE-1.X/licenses/


    - SINGULARITY_EXEC:

    	Path of the singularity container that contains the pele executable. It can be configured in two ways:

    	  1. Can be specified by input.yaml flag

    	  .. code-block:: yaml

    	  	singularity_exec: /path/to/singularity/container/


    	  2. Can be specified by environment variable:

    	  .. code-block:: bash

    		export PELE_LICENSE=/path/to/singularity/container/