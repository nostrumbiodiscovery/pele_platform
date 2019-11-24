import subprocess
import pele_platform.Utilities.Parameters.pele_env as pele
import os
import glob
import AdaptivePELE.analysis.splitTrajectory as st
import pele_platform.constants as cs
import pele_platform.Utilities.Helpers.template_builder as tb
import pele_platform.Adaptive.simulation as sim

def run(args):
    env = sim.run_adaptive(args) 
    if args.restart in ["all", "adaptive", "score"]:
        env = pele.EnviroBuilder.build_env(args)
        split_traj_dirs = split_trajectory(env)
        rescore(split_traj_dirs, env)


def split_trajectory(paths):
    files = glob.glob(os.path.join(paths.adap_ex_output, "*/traj*.*"))
    epoch_files = [report for report in files if(os.path.basename(os.path.dirname(report)).isdigit())]
    for file in epoch_files:
        output_dir = os.path.join(paths.pele_dir, "ini_str", os.path.splitext(os.path.basename(file))[0])
        st.main(output_dir, [file, ], paths.topology, None, template=None)
        yield output_dir


def rescore(split_trajs_dirs, env):
    for split_traj_dir in split_trajs_dirs:
        snapshots = glob.glob(os.path.join(split_traj_dir, "*"))
        symbolic_links_to_structs(snapshots, env)
    score_glide(env)


def symbolic_links_to_structs(snapshots, env):
     for snapshot in snapshots:
        try:
            os.symlink(os.path.abspath(snapshot), os.path.join(env.glide_structs, os.path.basename(snapshot)))
        except OSError:
            pass


def score_glide(env):
    schrodinger_exec = os.path.join(cs.SCHRODINGER, "run")
    tb.TemplateBuilder(env.glide_template, {"INPUT": env.glide_structs, "PRECISION": env.precision_glide})
    command = "{0} xglide.py {1} -HOST {2}:{3}".format(schrodinger_exec, env.glide_template, "Calculon_slurm", env.cpus)
    print(command)
    #subprocess.call(command.split())
