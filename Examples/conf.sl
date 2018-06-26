#!/bin/bash
#SBATCH -J PELE_MPI
#SBATCH --output=mpi_%j.out
#SBATCH --error=mpi_%j.err
#SBATCH --qos=debug
#SBATCH --ntasks=5

module load intel mkl impi python/2.7.13

export PYTHONPATH=/gpfs/projects/bsc72/MSM_XTC/:/gpfs/projects/bsc72/lib/site-packages/:$PYTHONPATH

#python -m MSM_PELE.main  PR_1A28_xray_-_minimized.pdb STR Z --cpus 5 --test -wf /home/bsc72/bsc72893/test_MSM/STR_Pele_9
python -m MSM_PELE.main  PR_1A28_xray_-_minimized.pdb STR Z --cpus 5 --test 
