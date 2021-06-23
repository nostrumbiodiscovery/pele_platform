#!/bin/bash
#SBATCH --job-name=PELE
#SBATCH --output=PELE.out
#SBATCH --error=PELE.err
#SBATCH --ntasks=5    # change this value to match CPUs in input.yaml
#SBATCH --qos=debug
#SBATCH --time=00-01:00:00

module purge
export PELE="/gpfs/projects/bsc72/PELE++/mniv/V1.6.2-b1/"
export SCHRODINGER="/gpfs/projects/bsc72/SCHRODINGER_ACADEMIC"
export PATH=/gpfs/projects/bsc72/conda_envs/platform/1.5.1/bin:$PATH
module load intel mkl impi gcc # 2> /dev/null
module load impi
module load boost/1.64.0
/gpfs/projects/bsc72/conda_envs/platform/1.5.1/bin/python3.8  -m pele_platform.main input.yaml