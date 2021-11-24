=======================
Input files preparation
=======================

PELE Platform essentially needs a minimum of two files in order to run:
    - ``system.pdb``: a PDB file that contains the structure to simulate with PELE.
      It typically contains a protein and a small molecule.
    - ``input.yaml``: a yaml file containing a collection of parameters to
      set up the simulation.


Below we list several instructions to correctly prepare input PDB files for PELE:

.. toctree::
   :maxdepth: 2

   system_preparation.rst


On the other hand, the list of all parameters that can be specified in the
``input.yaml`` file can be found below:

.. toctree::
   :maxdepth: 2

   yaml.rst
