#!/bin/bash
#SBATCH -J PELE_MPI
#SBATCH --output=mpi_%j.out
#SBATCH --error=mpi_%j.err
#SBATCH --ntasks=11
#SBATCH --mem-per-cpu=1000

module purge

#Try 2015 GNU

#module load Python
#module load PELE/1.5.0.2524

#Try 2018 GNU

unset PYTHONPATH
unset LD_LIBRARY_PATH

module load Python/2.7.14-foss-2018a  Boost/1.66.0-foss-2018a JsonCpp/1.8.4-foss-2018a Crypto++/6.1.0-intel-2018a patchelf/0.9-foss-2018a wjelement/1.3-foss-2018a CMake/3.7.2-GCCcore-6.3.0 GCCcore/6.4.0 OpenMPI
module load Python/2.7.14-foss-2018a
export PYTHONPATH=/sNow/easybuild/centos/7.4.1708/Skylake/software/PyMOL/2.2.0_0/lib/python2.7/site-packages/:$PYTHONPATH
export PYTHONPATH=/work/NBD_Utilities/PELE/PELE_Softwares/:/work/NBD_Utilities/PELE/PELE_Softwares/PELE_Platform/:$PYTHONPATH

#KINASE TEST
#python -m PELE_Platform.main /work/NBD_Utilities/PELE/PELE_Softwares/PELE_Platform/Examples/Kinase/1_3ZON_complex.pdb IK1 Z --test --hbond A:690:_H__ Z:1:_O2_

#MSM TEST PDB
#python /work/NBD_Utilities/PELE/PELE_Softwares/PELE_Platform/main.py /work/NBD_Utilities/PELE/PELE_Softwares/PELE_Platform/Examples/Msm/PR_1A28_xray_-_minimized.pdb STR Z --test --msm 

#MSM TEST MAE
#python /work/NBD_Utilities/PELE/PELE_Softwares/PELE_Platform/main.py /work/NBD_Utilities/PELE_Softwares/PELE_Platform/Examples/Msm/PR_1A28_protein.pdb AS4 Z --mae_lig AS4.mae --test --msm --restart adaptive

#ADAPTIVE TEST 
#python /work/NBD_Utilities/PELE/PELE_Softwares/PELE_Platform/main.py /work/NBD_Utilities/PELE/PELE_Softwares/PELE_Platform/Examples/Msm/PR_1A28_xray_-_minimized.pdb STR Z --adaptive /work/NBD_Utilities/PELE/PELE_Softwares/PELE_Platform/Examples/Adaptive/adaptive.conf --pele /work/NBD_Utilities/PELE/PELE_Softwares/PELE_Platform/Examples/Adaptive/pele_file_4Adaptive.conf --cpus 2 --pdb

#ADAPTIVE TEST EXTERNAL TEMPLATES
#python /work/NBD_Utilities/PELE_Softwares/PELE_Platform/main.py /work/NBD_Utilities/PELE_Softwares/PELE_Platform/Examples/Msm/PR_1A28_xray_-_minimized.pdb STR Z --adaptive /work/NBD_Utilities/PELE_Softwares/PELE_Platform/Examples/Adaptive/adaptive.conf --pele /work/NBD_Utilities/PELE_Softwares/PELE_Platform/Examples/Adaptive/pele_file_4Adaptive.conf --cpus 2 --pdb --template /work/NBD_Utilities/PELE_Softwares/PELE_Platform/Examples/Adaptive/strz --rotamers /work/NBD_Utilities/PELE_Softwares/PELE_Platform/Examples/Adaptive/STR.rot.assign

#ADAPTIVE OUT IN TEST
#python -m PELE_Platform.main /work/NBD_Utilities/PELE/PELE_Softwares/PELE_Platform/Examples/Kinase/1_3ZON_complex.pdb IK1 Z --test --out_in

#ADAPTIVE GLOBAL TEST
#python -m PELE_Platform.main /work/NBD_Utilities/PELE/PELE_Softwares/PELE_Platform/Examples/Kinase/1_3ZON_complex.pdb IK1 Z --full --cpus 11

#ADAPTIVE INDUCE FIT TEST
#python -m PELE_Platform.main /work/NBD_Utilities/PELE/PELE_Softwares/PELE_Platform/Examples/Msm/PR_1A28_xray_-_minimized.pdb STR Z --test --induce_fit

#FRAG PELE TEST NOOO
#python -m PELE_Plaform.grow_for_pele -cp /work/NBD_Utilities/PELE/PELE_Softwares/PELE_Platform/Examples/Frag/181l_prepared.pdb -fp /work/NBD_Utilities/PELE/PELE_Softwares/PELE_Platform/Examples/Frag/methyl.pdb -ca C4 -fa C1

#ADAPTIVE IN OUT FORCED TEST
#python -m PELE_Platform.main /work/NBD_Utilities/PELE/PELE_Softwares/PELE_Platform/Examples/Msm/PR_1A28_xray_-_minimized.pdb STR Z --in_out --test

#ADAPTIVE IN OUT SOFT TEST NOOO
export PYTHONPATH=/home/dsoler/adaptive_types/md/AdaptivePELE/:$PYTHONPATH
#python -m PELE_Platform.main /work/NBD_Utilities/PELE/PELE_Softwares/PELE_Platform/Examples/Msm/PR_1A28_xray_-_minimized.pdb STR Z --in_out_soft --test

#ADAPTIVE WATER TEST
module purge
module load icc/2018.1.163-GCC-6.4.0-2.28 intel/2018a imkl/2018.1.163-iimpi-2018a  impi/2018.1.163-iccifort-2018.1.163-GCC-6.4.0-2.28 Python/2.7.14-intel-2018a  Boost/1.66.0-intel-2018a  patchelf/0.9-intel-2018a wjelement/1.3-foss-2018a JsonCpp/1.8.4-intel-2018a CMake/3.9.5-GCCcore-6.4.0 GCC/6.4.0-2.28 impi/2018.1.163-iccifort-2018.1.163-GCC-6.4.0-2.28
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so
#python -m PELE_Platform.main /work/NBD_Utilities/PELE/PELE_Softwares/PELE_Platform/Examples/water/water_processed.pdb HOH M --water_exp M:1 --cpus 5

#ADAPTIVE WATER LIGAND TEST
#python -m PELE_Platform.main /work/NBD_Utilities/PELE/PELE_Softwares/PELE_Platform/Examples/water/hit1_complex_processed.pdb LIG L --water_lig M:1 --cpus 5

