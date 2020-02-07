#!/bin/bash
#BSUB -J PELEPlatform
#BSUB -W 00:50
#BSUB -q bsc_debug
#BSUB -eo output_%J.err
#BSUB -oo output_%J.out
#BSUB -n 11

module purge
module load intel gcc openmpi/1.8.1 boost/1.63.0 python/2.7.3 MKL/11.3 GTK+3/3.2.4

export PYTHONPATH=/gpfs/projects/bsc72/PELEPlatform/1.2.3/:/gpfs/projects/bsc72/adaptiveSampling/bin_nord/v1.6.2/::/gpfs/projects/bsc72/PELEPlatform/external_deps/:$PYTHONPATH
export PYTHONPATH=/gpfs/projects/bsc72/lib/site-packages_mn3:$PYTHONPATH

export MPLBACKEND=Agg
export OMPI_MCA_coll_hcoll_enable=0
export OMPI_MCA_mtl=^mxm

python -m pele_platform.main input.yaml

