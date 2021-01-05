############################################################
#                                                          #
#   RUNNING PRODUCTION ON THE OFFICE MACHINE (nbdcalc01)   #
#                                                          #
############################################################

# export when running tests from nbdcalc01 (it cannot see cluster files)
export LC_ALL=C; unset LANGUAGE
export PELE=/scratch/PELE-repo/
export PATH="/usr/lib64/openmpi/bin/":$PATH
export SCHRODINGER=/opt/schrodinger2020-1/
export LD_LIBRARY_PATH=/scratch/PELE-compilation/PELE-dependencies/boost-1.52.0/lib/

# easy way to remove files from previous test runs
rm -r LIG_Pele* SB4_Pele* STR_Pele* AS4_Pele* IK1_Pele* allosteric NOR_solvent_OBC/ PCA_result/ *.pdb *.log API_Pele* 1w7h_preparation_structure_2w_processed*
#production run
/shared/work/NBD_Utilities/PELE/PELE_Softwares/PelePlatform/depend/bin/python -m pele_platform.main input.yaml

