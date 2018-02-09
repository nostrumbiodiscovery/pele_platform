import socket
machine = socket.getfqdn()

if "bsc.mn" in machine:
    SCHRODINGER = "/gpfs/projects/bsc72/SCHRODINGER_ACADEMIC"
    PELE = "/gpfs/projects/bsc72/PELE++/mniv/rev12455"
    ADAPTIVE = "/gpfs/projects/bsc72/adaptiveSampling/bin/v1.4.1"
    MPIRUN = "/apps/INTEL/2017.4/impi/2017.3.196/bin64"
else:
    SCHRODINGER = "/sNow/easybuild/centos/7.4.1708/Skylake/software/schrodinger2017-4/"
    PELE = "/sNow/easybuild/centos/7.4.1708/Skylake/software/PELE/1.5.0.2524/"
    ADAPTIVE = "/home/dsoler/repos/AdaptivePELE/"
    MPIRUN = "/sNow/easybuild/centos/7.4.1708/Skylake/software/OpenMPI/1.8.4-GCC-4.9.2/bin"

