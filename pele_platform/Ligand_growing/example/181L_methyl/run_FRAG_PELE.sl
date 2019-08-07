#!/bin/bash
#SBATCH -J FRAG_PELE_MPI
#SBATCH --output=MSM_%j.out
#SBATCH --error=MSM_%j.err
#SBATCH --ntasks=60
#SBATCH --mem-per-cpu=1G

OMP_NUM_THREADS=1

module purge
module load GCCcore/6.3.0
module load Boost/1.66.0-foss-2018a JsonCpp/1.8.4-foss-2018a Crypto++/6.1.0-intel-2018a patchelf/0.9-foss-2018a wjelement/1.3-foss-2018a CMake/3.7.2-GCCcore-6.3.0 GCCcore/6.4.0 OpenMPI
module load Python/3.6.4-foss-2018a
export PYTHONPATH=/home/dsoler/repos/AdaptivePELE/:$PYTHONPAT
python ../../grow_for_pele.py -cp 181l_prepared.pdb -fp methyl.pdb -ca C4 -fa C1

