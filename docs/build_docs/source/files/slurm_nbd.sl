#!/bin/bash
#SBATCH -J PELE_MPI
#SBATCH --output=mpi_%j.out
#SBATCH --error=mpi_%j.err
#SBATCH --ntasks=60        # Change this value to match CPUs in input.yaml
#SBATCH --mem-per-cpu=2000

##################################################################
module purge
export SCHRODINGER="/sNow/easybuild/centos/7.4.1708/Skylake/software/schrodinger2017-4/"
export PELE="/work/NBD_Utilities/PELE/PELE_Softwares/bin/PELE1.6/"
export LC_ALL=C; unset LANGUAGE
unset PYTHONPATH
module load impi/2018.1.163-iccifort-2018.1.163-GCC-6.4.0-2.28 Boost/1.66.0-intel-2018a wjelement/1.3-intel-2018a
module load Crypto++/6.1.0-intel-2018a OpenBLAS/0.2.20-GCC-6.4.0-2.28
module load Python/3.7.0-foss-2018a
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so
export PYTHONPATH=/work/NBD_Utilities/PELE/PELE_Softwares/PelePlatform/pele_platform/:/work/NBD_Utilities/PELE/PELE_Softwares/PelePlatform/pele_platform_dependencies/:$PYTHONPATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/dsoler/pele_deps/boost_1_52/lib/
module load RDKit/2019.09.3-foss-2019b-Python-3.7.4
##################################################################

python -m pele_platform.main input.yaml