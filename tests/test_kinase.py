import pytest
import os
import pele_platform.constants.constants as cs
import pele_platform.main as main

test_path = os.path.join(cs.DIR, "Examples")





KINASE_ARGS = [os.path.join(test_path, "Kinase/1_3ZON_complex.pdb"), "IK1", "Z", "--test", "--hbond", "A:690:_H__", "Z:1:_O2_"]

@pytest.mark.parametrize("ext_args", [
                         (KINASE_ARGS),
                         ])
def test_kinases(ext_args):
    arguments = main.parseargs(ext_args)
    main.set_software_to_use(arguments)
    main.Launcher(arguments).launch()



BIAS_ARGS = [os.path.join(test_path, "Bias/complex.pdb"), "LIG", "L", "--test", "--bias", "--atom_dist", "L:275:N", "C:89:CA", 
    "--pele_steps", "1", "--iterations", "1"] 
@pytest.mark.parametrize("ext_args", [
                         (BIAS_ARGS),
                         ])
def test_bias(ext_args):
    arguments = main.parseargs(ext_args)
    main.set_software_to_use(arguments)
    main.Launcher(arguments).launch()


EXTERNAL_CONFILES_ARGS= [os.path.join(test_path, "Msm/PR_1A28_xray_-_minimized.pdb"), "STR", "Z", "--cpus", "2", "--pele", os.path.join(test_path, "Adaptive/pele_file_4Adaptive.conf"), "--adaptive", os.path.join(test_path, "Adaptive/adaptive.conf"), "--template", os.path.join(test_path, "Adaptive/strz"), "--rotamers", os.path.join(test_path, "Adaptive/STR.rot.assign")] 
@pytest.mark.parametrize("ext_args", [
                         (EXTERNAL_CONFILES_ARGS),
                         ])
def test_external_confiles(ext_args):
    arguments = main.parseargs(ext_args)
    main.set_software_to_use(arguments)
    main.Launcher(arguments).launch()


#python -m PELE_Platform.main /work/NBD_Utilities/PELE/PELE_Softwares/PELE_Platform/Examples/Kinase/1_3ZON_complex.pdb IK1 Z --test --out_in
OUT_IN_ARGS = [os.path.join(test_path, "Kinase/1_3ZON_complex.pdb"), "IK1", "Z", "--bias", "--out_in", "--solvent", "OBC", "--iterations", "1",
"--pele_steps", "1", "--cpus", "3", "--report_name", "report_sim1", "--traj_name", "traj_sim1"]
@pytest.mark.parametrize("ext_args", [
                         (OUT_IN_ARGS),
                         ])
def test_out_in(ext_args):
    arguments = main.parseargs(ext_args)
    main.set_software_to_use(arguments)
    main.Launcher(arguments).launch()

#python -m PELE_Platform.main /work/NBD_Utilities/PELE/PELE_Softwares/PELE_Platform/Examples/Kinase/1_3ZON_complex.pdb IK1 Z --full --cpus 11
GLOBAL_ARGS = [os.path.join(test_path, "Kinase/1_3ZON_complex.pdb"), "IK1", "Z", "--full", "--cpus", "11"]
@pytest.mark.parametrize("ext_args", [
                         (GLOBAL_ARGS),
                         ])
def test_global(ext_args):
    arguments = main.parseargs(ext_args)
    main.set_software_to_use(arguments)
    main.Launcher(arguments).launch()

#python -m PELE_Platform.main /work/NBD_Utilities/PELE/PELE_Softwares/PELE_Platform/Examples/Msm/PR_1A28_xray_-_minimized.pdb STR Z --test --induce_fit --atom_dist Z:1:H8 A:1:CA 
INDUCE_FIT_ARGS = [os.path.join(test_path, "Msm/PR_1A28_xray_-_minimized.pdb"), "STR", "Z", "--induce_fit", "--test"]
@pytest.mark.parametrize("ext_args", [
                         (INDUCE_FIT_ARGS),
                         ])
def test_induce_fit(ext_args):
    arguments = main.parseargs(ext_args)
    main.set_software_to_use(arguments)
    main.Launcher(arguments).launch()

#python -m PELE_Platform.main /work/NBD_Utilities/PELE/PELE_Softwares/PELE_Platform/Examples/Msm/PR_1A28_xray_-_minimized.pdb STR Z --in_out --test
EXIT_ARGS = [os.path.join(test_path, "Msm/PR_1A28_xray_-_minimized.pdb"), "STR", "Z", "--in_out", "--test", "--epsilon", "0.19", 
    "--pele_steps", "1", "--iterations", "1"] 
@pytest.mark.parametrize("ext_args", [
                         (EXIT_ARGS),
                         ])
def test_exit(ext_args):
    arguments = main.parseargs(ext_args)
    main.set_software_to_use(arguments)
    main.Launcher(arguments).launch()

#python -m PELE_Platform.main /work/NBD_Utilities/PELE/PELE_Softwares/PELE_Platform/Examples/water/water_processed.pdb HOH M --water_exp M:1 --cpus 5
WATER_ARGS = [os.path.join(test_path, "water/water_processed.pdb"), "HOH", "M", "--water_exp", "M:1", "--cpus", "5", "--pele_steps", "1", "--iterations", "1"]

@pytest.mark.parametrize("ext_args", [
                         (WATER_ARGS),
                         ])
def test_water(ext_args):
    arguments = main.parseargs(ext_args)
    main.set_software_to_use(arguments)
    main.Launcher(arguments).launch()

#python -m PELE_Platform.main /work/NBD_Utilities/PELE/PELE_Softwares/PELE_Platform/Examples/water/hit1_complex_processed.pdb LIG L --water_lig M:1 --cpus 5 --water_center 14.505 -25.051  -1.524

LIG_WATER_ARGS = [os.path.join(test_path, "water/hit1_complex_processed.pdb"), "LIG", "L", "--water_lig", "M:1", "--cpus", "5", "--water_center", "14.505", "-25.051", "-1.524", "--pele_steps", "1", "--iterations", "1"]

@pytest.mark.parametrize("ext_args", [
                         (LIG_WATER_ARGS),
                         ])
def test_lig_water(ext_args):
    arguments = main.parseargs(ext_args)
    main.set_software_to_use(arguments)
    main.Launcher(arguments).launch()


#python /work/NBD_Utilities/PELE/PELE_Softwares/PELE_Platform/main.py /work/NBD_Utilities/PELE/PELE_Softwares/PELE_Platform/Examples/Msm/PR_1A28_xray_-_minimized.pdb STR Z --test --msm 
MSM_PDB_ARGS = [os.path.join(test_path, "Msm/PR_1A28_xray_-_minimized.pdb"), "STR", "Z", "--test", "--msm"]

@pytest.mark.parametrize("ext_args", [
                         (MSM_PDB_ARGS),
                         ])
def test_msm_pdb(ext_args):
    arguments = main.parseargs(ext_args)
    main.set_software_to_use(arguments)
    main.Launcher(arguments).launch()


#python /work/NBD_Utilities/PELE/PELE_Softwares/PELE_Platform/main.py /work/NBD_Utilities/PELE_Softwares/PELE_Platform/Examples/Msm/PR_1A28_protein.pdb AS4 Z --mae_lig AS4.mae --test --msm --restart adaptive

MSM_MAE_ARGS = [os.path.join(test_path, "Msm/PR_1A28_protein.pdb"), "AS4", "Z", "--mae_lig", os.path.join(test_path, "Msm/AS4_INIT.mae"), "--test", "--msm"]

@pytest.mark.parametrize("ext_args", [
                         (MSM_MAE_ARGS),
                         ])
def test_msm_mae(ext_args):
    arguments = main.parseargs(ext_args)
    main.set_software_to_use(arguments)
    main.Launcher(arguments).launch()


PREPWIZARD_ARGS = [os.path.join(test_path, "Bias/complex.pdb"), "LIG", "L", "--test", "--bias", "--atom_dist", "L:275:N", "C:89:CA", 
    "--pele_steps", "1", "--iterations", "1", "--prepwizard"] 
@pytest.mark.parametrize("ext_args", [
                         (PREPWIZARD_ARGS),
                         ])
def test_prepwizard(ext_args):
    arguments = main.parseargs(ext_args)
    main.set_software_to_use(arguments)
    main.Launcher(arguments).launch()

if __name__ == "__main__":
    #test_kinases(KINASE_ARGS)
    #test_bias(BIAS_ARGS)
    #test_external_confiles(EXTERNAL_CONFILES_ARGS)
    #test_out_in(OUT_IN_ARGS)
    #test_global(GLOBAL_ARGS)
    #test_induce_fit(INDUCE_FIT_ARGS)
    #test_exit(EXIT_ARGS)
    #NOT WORKING WITH NEW
    #test_water(WATER_ARGS)
    #NOT WORKING WITH NEW
    #test_lig_water(LIG_WATER_ARGS)
    test_prepwizard(PREPWIZARD_ARGS)
    #test_msm_pdb(MSM_PDB_ARGS)
    #test_msm_mae(MSM_MAE_ARGS)
