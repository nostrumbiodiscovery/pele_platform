import os
import argparse
import schrodinger.structure as st

# Here because otherwise need to charge helpers with PPP and schrodinger does not have that module


class cd:
    """Context manager for changing the current working directory"""

    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)


def pdb_to_mae(fname, schr_path, mae_output_file=None, remove=False):
    directory = os.path.dirname(fname)
    fname = os.path.basename(fname)

    with cd(directory):
        file_info = os.path.splitext(fname)[0].split("_")
        try:
            epoch, trajectory, snapshot = [int(a) for a in file_info[0].split(".")]  # for BestStructs
        except:
            epoch, trajectory, snapshot = [int(a) for a in file_info[1].split(".")]  # clusters have different filenames
        binding_energy, = [float(a.replace("BindEner", "")) for a in file_info if "BindEner" in a]

        properties = {
            "BindingEnergy": binding_energy,
            "trajectory": trajectory,
            "snapshot": snapshot,
            "epoch": epoch,
        }

        traj = "{}_{}_{}".format(epoch, trajectory, snapshot)
        pele_energy = binding_energy
        cmd = (
            "$SCHRODINGER/utilities/prepwizard -rehtreat -noepik -disulfides -noimpref -nometaltreat -noprotassign -nopropka -WAIT %s %s.mae"
            % (fname, traj)
        )
        print(cmd)
        os.system(cmd)
        struct = next(st.StructureReader("{}.mae".format(traj)))
        struct.property["r_user_PELE_energy"] = properties["BindingEnergy"]
        struct.property["r_user_PELE_epoch"] = properties["epoch"]
        struct.property["r_user_PELE_traj"] = properties["trajectory"]
        struct.property["r_user_PELE_snapshot"] = properties["snapshot"]
        struct.title = "%s_BindEn_%.2f" % (traj, pele_energy)
        struct.write("{}.mae".format(traj))


def add_args(parser):
    parser.add_argument("inputfile", type=str, help="Pdb input file")
    parser.add_argument("--schr", type=str, help="Schrodinger root path")
    parser.add_argument(
        "--remove", action="store_true", help="Remove input file at exit"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Write MAE and SDF files with certain properties"
    )
    add_args(parser)
    args = parser.parse_args()
    pdb_to_mae(args.inputfile, args.schr, remove=args.remove)
