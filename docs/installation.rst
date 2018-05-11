============
Installation
============

.. toctree::
   :maxdepth: 2

MSM PELE installation requires several dependencies they must be set
before practice. Furthermore MSM PELE requires Schrodinger and PELE installation.

Conda Installation (Recommended)
---------------------------------

**Create conda**::

  conda create -y -n py27 python=2.7 anaconda

**Activate Env**::

  source activate py27

**Install dependencies**::

  ~/conda/envs/py27/bin/pip install msmtools

  ~/conda/envs/py27/bin/pip install pyemma

  ~/conda/envs/py27/bin/pip install prody==1.8.2

  ~/conda/envs/py27/bin/pip install sklearn

 
MSM PELE Configuration
-----------------------

**Enviromental variables**::
 
  export PYTHONPATH=$PYTHONPATH:'/path/to/MSM_PELE/'

**Changing local path**::

  Change the public variables path under **MSM_PELE/constants.py** to the ones local to your machine

**Pyemma config**::

  $ python

  > import pyemma

  > pyemma.config.used_filenames

You will recieve a list with files where it is possible to find configuration from pyemma. You need to go over these files and change::

    show_progress_bars = Flase

    mute= True
