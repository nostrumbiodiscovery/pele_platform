.. pele platform master file, created by
   sphinx-quickstart on Mon Dec  4 11:58:09 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Pele Platform
===========================================

Pele platform is a python module to automatically
launch PELE and AdaptivePELE. We will support  MSMPELE and FragPELE in next releases. 
It is built as a python layer on top of all PELE functionalities that can be used as
a backend for a potential future GUI.

Github : https://github.com/NostrumBioDiscovery/pele_platform

Requirements
===================

- Academic Schordinger

- Pymol Python

  .. code-block:: bash

    yum install gcc gcc-c++ kernel-devel python-devel tkinter python-pmw glew-devel freeglut-devel libpng-devel freetype-devel libxml2-devel glm-devel

    git clone https://github.com/schrodinger/pymol-open-source.git

    cd pymol-open-source

    python setup.py install (with same python you will use later)


Documentation
===================

.. toctree::
   installation/index.rst

.. toctree::
   cheatsheet/index.rst

.. toctree::
   documentation/index.rst

.. toctree::
   changelog/index.rst

