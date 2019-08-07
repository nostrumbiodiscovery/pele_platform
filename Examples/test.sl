#!/bin/bash
#SBATCH -J PELE_MPI
#SBATCH --output=mpi_%j.out
#SBATCH --error=mpi_%j.err
#SBATCH --ntasks=11
#SBATCH --mem-per-cpu=1000

module purge
unset PYTHONPATH
unset LD_LIBRARY_PATH
module load impi 
module load Python/2.7.14-foss-2018a
export PYTHONPATH=/sNow/easybuild/centos/7.4.1708/Skylake/software/PyMOL/2.2.0_0/lib/python2.7/site-packages/:$PYTHONPATH
export PYTHONPATH=~/development/pele_platform/pele_platform:~/development/pele_platform:$PYTHONPATH

#KINASE TEST
#python -m pele_platform.main /work/NBD_Utilities/PELE/PELE_Softwares/PELE_Platform/Examples/Kinase/1_3ZON_complex.pdb IK1 Z --test --hbond A:690:_H__ Z:1:_O2_

#MSM TEST PDB
#python /work/NBD_Utilities/PELE/PELE_Softwares/pele_platform/main.py /work/NBD_Utilities/PELE/PELE_Softwares/pele_platform/Examples/Msm/PR_1A28_xray_-_minimized.pdb STR Z --test --msm 

#MSM TEST MAE
#python /work/NBD_Utilities/PELE/PELE_Softwares/pele_platform/main.py /work/NBD_Utilities/PELE_Softwares/pele_platform/Examples/Msm/PR_1A28_protein.pdb AS4 Z --mae_lig AS4.mae --test --msm --restart adaptive

#YOUR OWN ADAPTIVE PELE
#python -m pele_platform.main /work/NBD_Utilities/PELE/PELE_Softwares/PELE_Platform/Examples/Msm/PR_1A28_xray_-_minimized.pdb STR Z --adaptive /work/NBD_Utilities/PELE/PELE_Softwares/PELE_Platform/Examples/Adaptive/adaptive.conf --pele /work/NBD_Utilities/PELE/PELE_Softwares/PELE_Platform/Examples/Adaptive/pele_file_4Adaptive.conf --cpus 2 --pdb --template /work/NBD_Utilities/PELE/PELE_Softwares/PELE_Platform/Examples/Adaptive/strz --rotamers /work/NBD_Utilities/PELE/PELE_Softwares/PELE_Platform/Examples/Adaptive/STR.rot.assign

#ADAPTIVE OUT IN TEST
#python -m pele_platform.main /work/NBD_Utilities/PELE/PELE_Softwares/PELE_Platform/Examples/Kinase/1_3ZON_complex.pdb IK1 Z --test --out_in

#ADAPTIVE GLOBAL TEST
#python -m pele_platform.main /work/NBD_Utilities/PELE/PELE_Softwares/PELE_Platform/Examples/Kinase/1_3ZON_complex.pdb IK1 Z --full --cpus 11

#ADAPTIVE INDUCE FIT TEST
#python -m pele_platform.main /work/NBD_Utilities/PELE/PELE_Softwares/PELE_Platform/Examples/Msm/PR_1A28_xray_-_minimized.pdb STR Z --test --induce_fit --atom_dist Z:1:H8 A:1:CA 

#ADAPTIVE IN OUT 
#python -m pele_platform.main /work/NBD_Utilities/PELE/PELE_Softwares/PELE_Platform/Examples/Msm/PR_1A28_xray_-_minimized.pdb STR Z --in_out --cpus 11

#ADAPTIVE WATER TEST
#python -m pele_platform.main /work/NBD_Utilities/PELE/PELE_Softwares/PELE_Platform/Examples/water/water_processed.pdb HOH M --water_exp M:1 --cpus 5

#ADAPTIVE WATER LIGAND TEST NOOO
#python -m pele_platform.main /work/NBD_Utilities/PELE/PELE_Softwares/PELE_Platform/Examples/water/hit1_complex_processed.pdb LIG L --water_lig M:1 --cpus 5 --water_center 14.505 -25.051  -1.524

