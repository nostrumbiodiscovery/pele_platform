#!/bin/bash
#SBATCH -J PELE_MPI
#SBATCH --time=02:00:00
#SBATCH --ntasks=3
#SBATCH --mem-per-cpu=1G

module load Cython PELE/1.5.0.2524

python -m AdaptivePELE.adaptiveSampling $CONFILE


