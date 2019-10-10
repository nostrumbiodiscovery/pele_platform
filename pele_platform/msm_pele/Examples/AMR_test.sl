#!/bin/bash
#SBATCH -J PELE_MPI
#SBATCH --output=mpi_%j.out
#SBATCH --error=mpi_%j.err
#SBATCH --ntasks=10
#SBATCH --qos=debug
module purge
module load intel mkl impi python/2.7.13 boost/1.64.0_py2 gcc 
python -m MSM_PELE.main URO_INIT.pdb AMR B --mae_lig 1F5L.mae --test  --iterations 3
