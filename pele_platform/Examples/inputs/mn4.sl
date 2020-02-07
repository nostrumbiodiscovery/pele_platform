#!/bin/bash

#SBATCH --job-name=PELEPlatformTests
#SBATCH --output=output_%j.out
#SBATCH --error=output_%j.err
#SBATCH --ntasks=10
#SBATCH --time=00-02:00:00
#SBATCH -D .
#SBATCH --qos=debug

#############################NO CHANGE###########################
module purge #2> /dev/null
module load intel mkl impi gcc # 2> /dev/null
module load python/2.7.13 # 2> /dev/null
module load boost/1.64.0_py2 # 2> /dev/null
module load python/2.7.13 glew/2.1.0 glm/0.9.9.6 qt/5.12.1
export PYTHONPATH="/gpfs/projects/bsc72/PELEPlatform/1.2.3:/gpfs/projects/bsc72/adaptiveSampling/bin/v1.6.2/:/gpfs/projects/bsc72/lib_msm/site-packages/:/gpfs/projects/bsc72/lib/site-packages:/gpfs/projects/bsc72/PELEPlatform/external_deps/"
############################CHANGE##########################

python -m pele_platform.main input.yaml
