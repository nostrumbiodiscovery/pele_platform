#!/bin/bash
#SBATCH -J PELE_MPI
#SBATCH --output=mpi_%j.out
#SBATCH --error=mpi_%j.err
#SBATCH --ntasks=11
#SBATCH --mem-per-cpu=1000

#############################NO CHANGE###########################
module purge
unset PYTHONPATH
unset LD_LIBRARY_PATH
module load impi 
module load Python/2.7.14-foss-2018a
export PYTHONPATH=/sNow/easybuild/centos/7.4.1708/Skylake/software/PyMOL/2.2.0_0/lib/python2.7/site-packages/:$PYTHONPATH
export PYTHONPATH=/work/NBD_Utilities/PELE/PELE_Softwares/:/work/NBD_Utilities/PELE/PELE_Softwares/PELE_Platform/:$PYTHONPATH



############################CHANGE##########################


pytest /work/NBD_Utilities/PELE/PELE_Softwares/PELE_Platform/tests -vv > out.txt
