#!/bin/bash
#SBATCH -J PELE_MPI
#SBATCH --output=mpi_%j.out
#SBATCH --error=mpi_%j.err
############################CHANGE##########################
#SBATCH --ntasks=10                                      #Example--> #SBATCH --ntasks=250
############################CHANGE##########################
#SBATCH --mem-per-cpu=1000

#########USER VARIABLE#########

INPUT_NAME="input.yaml"


########IT TO CHANGE#####################

PELE_LICENSES="/path/to/dir/PELE/licenses/"
SCHRODINGER="/path/to/schordinger20XX/"
SCHR_LIC="/path/to/schodigner/license/dir/"
PATH_TO_PELE_IMG="/path/to/PELE/IMAGE/pele.img"

#Example:
#PELE_LICENSES="/work/NBD_Utilities/PELE/PELE_Softwares/bin/PELE1.6/licenses/"
#SCHRODINGER="/sNow/easybuild/centos/7.4.1708/Skylake/software/schrodinger2017-4/"
#SCHR_LIC="/opt/schrodinger/licenses/"
#PATH_TO_PELE_IMG="/home/dsoler/singularity/pele_platform/singularity/pele.img"

########################################

######NOT TO CHANGE#################

singularity exec -B "$PELE_LICENSES:/PELE-repo/licenses,$SCHRODINGER:/opt/,$SCHR_LIC:/schr_lic/" $PATH_TO_PELE_IMG  bash -c "source /etc/profile.d/01-locale-fix.sh; python -m pele_platform.main $INPUT_NAME"

####################################
