#!/bin/bash
#SBATCH -J PELE_MPI
#SBATCH --output=mpi_%j.out
#SBATCH --error=mpi_%j.err
#SBATCH --ntasks=5
#SBATCH --nodelist=node005
#SBATCH --mem-per-cpu=1000

module purge


#Try 2018 GNU
unset PYTHONPATH
unset LD_LIBRARY_PATH
module load GCCcore/6.3.0 
module load Python/2.7.14-foss-2018a  Boost/1.66.0-foss-2018a JsonCpp/1.8.4-foss-2018a Crypto++/6.1.0-intel-2018a patchelf/0.9-foss-2018a wjelement/1.3-foss-2018a CMake/3.7.2-GCCcore-6.3.0 GCCcore/6.4.0 OpenMPI
export PYTHONPATH=/home/dsoler/:/home/dsoler/PELE_Platform/:$PYTHONPATH

#KINASE TEST
python -m PELE_Platform.main Kinase/1_3ZON_complex.pdb IK1 Z --test --hbond A:690:_H__ Z:1:_O2_

#MSM TEST
#python -m PELE_Platform.main Msm/PR_1A28_xray_-_minimized.pdb STR Z --test --msm
#python -m PELE_Platform.main Msm/PR_1A28_protein.pdb AS4 Z --mae_lig AS4.mae --test --msm --restart adaptive

#ADAPTIVE TEST
#python -m PELE_Platform.main Adaptive/PR_1A28_xray_-_minimized.pdb STR Z --adaptive Adaptive/adaptive.conf --pele Adaptive/pele_file_4Adaptive.conf --cpus 2 --pdb



 
