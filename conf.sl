#!/bin/bash
#SBATCH -J PELE_MPI
#SBATCH --output=mpi_%j.out
#SBATCH --error=mpi_%j.err
#SBATCH --ntasks=70
#SBATCH --mem-per-cpu=1G

module load Python

#export LD_LIBRARY_PATH=<schrodinger_path>/mmshare-v4.0/lib/Linux-x86_64/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/sNow/easybuild/centos/7.4.1708/Skylake/software/schrodinger2017-4/mmshare-v4.0/lib/Linux-x86_64/:$LD_LIBRARY_PATH

#python /path/to/MSM_PELE/main.py <complex.pdb> <resname> <chain> --cpus X
python MSM_PELE/main.py PR_1A28_xray_-_minimized.pdb STR Z --cpus 5

