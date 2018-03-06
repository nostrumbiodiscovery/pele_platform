import socket
machine = socket.getfqdn()

if "bsc.mn" in machine:
    SCHRODINGER = "/gpfs/projects/bsc72/SCHRODINGER_ACADEMIC"
    PELE = "/gpfs/projects/bsc72/PELE++/mniv/rev12455"
    ADAPTIVE = "/gpfs/projects/bsc72/adaptiveSampling/bin/v1.4.2"
    MPIRUN = "/apps/INTEL/2017.4/impi/2017.3.196/bin64"
    LICENSE = "/gpfs/projects/bsc72/PELE++/license"
    MMSHARE = None

elif "mn.bsc" in machine:
    SCHRODINGER = "/gpfs/projects/bsc72/SCHRODINGER_ACADEMIC_NORD"
    MMSHARE = "/gpfs/projects/bsc72/SCHRODINGER_ACADEMIC_NORD/mmshare-v3.4/bin/Linux-x86_64"
    PELE = "/gpfs/projects/bsc72/PELE++/nord/rev12489"
    ADAPTIVE = "/gpfs/projects/bsc72/adaptiveSampling/bin_nord/v1.4.2"
    MPIRUN = "/apps/OPENMPI/1.8.1-mellanox/bin"
    LICENSE = "/gpfs/projects/bsc72/PELE++/license"

elif "bsccv" in machine:
    SCHRODINGER = "/data2/bsc72/SCHRODINGER_ACADEMIC"
    PELE = "/data/EAPM/PELE/PELE++/life/rev12489"
    ADAPTIVE = "/data2/bsc72/AdaptiveSampling/bin_nord/v1.4.2"
    MPIRUN = "/data2/apps/OPENMPI/1.6.1/bin"
    LICENSE = "/data/EAPM/PELE/PELE++/license"
    MMSHARE = "/data2/bsc72/SCHRODINGER_ACADEMIC/mmshare-v3.4/bin/Linux-x86_64"
    SCHRODINGER_PYTHON_LIBS = "/data2/bsc72/SCHRODINGER_ACADEMIC/mmshare-v3.4/lib/Linux-x86_64/lib/python2.7/site-packages/"
else:
    SCHRODINGER = "/sNow/easybuild/centos/7.4.1708/Skylake/software/schrodinger2017-4/"
    PELE = "/sNow/easybuild/centos/7.4.1708/Skylake/software/PELE/1.5.0.2524/"
    ADAPTIVE = "/home/dsoler/repos/AdaptivePELE/"
    MPIRUN = "/sNow/easybuild/centos/7.4.1708/Skylake/software/OpenMPI/1.8.4-GCC-4.9.2/bin"
    LICENSE = "/sNow/easybuild/centos/7.4.1708/Skylake/software/PELE/licenses/"
    MMSHARE = None
