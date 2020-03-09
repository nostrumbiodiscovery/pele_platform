#!/bin/bash
#SBATCH -J PELE_MPI
#SBATCH --output=mpi_%j.out
#SBATCH --error=mpi_%j.err
#SBATCH --ntasks=10
#SBATCH --mem-per-cpu=1000

#############################NO CHANGE###########################
module purge
unset PYTHONPATH
unset LD_LIBRARY_PATH
module load impi/2018.1.163-iccifort-2018.1.163-GCC-6.4.0-2.28 Boost/1.66.0-intel-2018a wjelement/1.3-intel-2018a
module load Crypto++/6.1.0-intel-2018a OpenBLAS/0.2.20-GCC-6.4.0-2.28
module load intel Python/2.7.14-intel-2018a
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so
export PYTHONPATH=/work/NBD_Utilities/PELE/PELE_Softwares/pele_platform_devel/:/work/NBD_Utilities/PELE/PELE_Softwares/adaptive_types/v1.6.2:$PYTHONPATH



############################CHANGE##########################

python -m pytest -s 
