import os
import glob
import MSM_PELE.msm as msm
import MSM_PELE.AdaptivePELE.analysis.splitTrajectory as st

def run(args):
    args.restart = "glide"
    env = msm.run(args)
    trajectory_snapshots = split_trajectory(env)
    #rescore_results = rescore(trajectory_snapshots)
    #plot(rescore_results)

def split_trajectory(paths):
    files = glob.glob(os.path.join(paths.adap_ex_output, "*/traj*.*"))
    epoch_files = [report for report in files if(os.path.basename(os.path.dirname(report)).isdigit())]
    print(epoch_files)
    for file in epoch_files:
        output_dir = os.path.join(paths.pele_dir, "ini_str", os.path.splitext(os.path.basename(file))[0])
        st.main(output_dir, files, paths.topology, None, template=None)

