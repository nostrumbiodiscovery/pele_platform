#!/bin/bash
#SBATCH -J PELE_MPI
#SBATCH --output=mpi_%j.out
#SBATCH --error=mpi_%j.err
#SBATCH --ntasks=90
#SBATCH --mem-per-cpu=1G


module purge

module load Python PELE

export PYTHONPATH=$PYTHONPATH:/path/to/MSM_PELE/

export LD_LIBRARY_PATH=/path/to/schrodinger2017-4/mmshare-v4.0/lib/Linux-x86_64/:$LD_LIBRARY_PATH

python /path/to/MSM_PELE/main.py  <PDB> <resname> <chain> --cpus <numcpus>
