#!/bin/bash
#SBATCH -J PELE_MPI_test
#SBATCH --output=mpi_%j.out
#SBATCH --error=mpi_%j.err
############################CHANGE##########################
#SBATCH --ntasks=50                                    #Example--> #SBATCH --ntasks=250
############################CHANGE##########################
#SBATCH --mem-per-cpu=1000

#############################NO CHANGE###########################
module purge
export SCHRODINGER="/sNow/easybuild/centos/7.4.1708/Skylake/software/schrodinger2017-4/"
export SCHRODINGER_PYTHONPATH="/sNow/easybuild/centos/7.4.1708/Skylake/software/schrodinger2017-4/internal/lib/python2.7/site-packages"
export PELE="/shared/work/NBD_Utilities/PELE/PELE_Softwares/bin/PELE1.6/"
export LC_ALL=C; unset LANGUAGE
unset PYTHONPATH
module load impi/2018.1.163-iccifort-2018.1.163-GCC-6.4.0-2.28 wjelement/1.3-intel-2018a
module load Crypto++/6.1.0-intel-2018a OpenBLAS/0.2.20-GCC-6.4.0-2.28
module load Python/3.7.0-foss-2018a
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so
export PYTHONPATH=/shared/work/NBD_Utilities/PELE/PELE_Softwares/PelePlatform/pele_platform/:/shared/work/NBD_Utilities/PELE/PELE_Softwares/PelePlatform/pele_platform_dependencies/:$PYTHONPATH
module load RDKit/2019.09.3-foss-2019b-Python-3.7.4
export LD_LIBRARY_PATH=/shared/work/NBD_Utilities/PELE/PELE_Softwares/local_deps/pele_deps/boost_1_52/lib:$LD_LIBRARY_PATH
#############################NO CHANGE###########################

export PYTHONPATH=/shared/home-nbdcalc01/dsoler/repos/pele_platform/:$PYTHONPATH

############################CHANGE##########################
#python -m pytest --cov=../ -s --cov-report=xml test*
python -m pele_platform.main input.yaml
############################CHANGE##########################
