import argparse
import os
import PELE_Platform.MSM.main as msm
import PELE_Platform.constants as cs
import PELE_Platform.Rescore.main as gl
import PELE_Platform.Rescore.simulation as ad


class Launcher():


    def __init__(self, args):
        self.args = args

    def launch(self):
        if args.software == "msm":
            if(args.clust > args.cpus and args.restart != "msm" and not args.test):
                raise ValueError(cs.CLUSTER_ERROR.format(args.cpus, args.clust))
            else:
                msm.run(args)

        elif args.software == "adaptive":
            ad.run_adaptive(args)

        elif args.software in ["glide", "out_in", "induce_fit"]:
            gl.run(args)


    @property
    def args(self):
        return self.args



def parseargs():
    parser = argparse.ArgumentParser(description='Run PELE Platform')
    parser.add_argument('system', type=str, help='complex to run pele on')
    parser.add_argument('residue', type=str, help='residue of the ligand to extract', default=cs.LIG_RES)
    parser.add_argument('chain', type=str, help='chain of the ligand to extract', default=cs.LIG_CHAIN)
    parser.add_argument("--mae_lig", type=str, help="ligand .mae file to include QM charges coming from jaguar")
    parser.add_argument("--box", type=str, help="Exploration box for Pele")
    parser.add_argument("--charge_ter", help="Charge protein terminals", action='store_true')
    parser.add_argument("--gaps_ter", help="Include TER when a possible gap is found", action='store_true')
    parser.add_argument("--clust", type=int, help="Numbers of clusters to start PELE's exploration with", default=cs.CLUSTERS)
    parser.add_argument('--forcefield', '-f', type=str, help='chain of the ligand to extract', default=cs.FORCEFIELD)
    parser.add_argument('--confile', type=str, help='your own pele configuration file', default=cs.PELE_CONFILE)
    parser.add_argument('--native', type=str, help='native file to compare RMSD to', default="")
    parser.add_argument('--cpus', type=int, help='number of processors', default=cs.CPUS)
    parser.add_argument("--core", type=int, help="Give one atom of the core section", default=-1)
    parser.add_argument("--mtor", type=int, help="Gives the maximum number of torsions allowed in each group.  Will freeze bonds to extend the core if necessary.", default=4)
    parser.add_argument("--n", type=int, help="Maximum Number of Entries in Rotamer File", default=1000)
    parser.add_argument("--clean", help="Whether to clean up all the intermediate files", action='store_true')
    parser.add_argument("--restart", type=str, help="Restart the platform from [all, pele, msm] with these keywords", default=cs.PLATFORM_RESTART)
    parser.add_argument("--gridres", type=str, help="Rotamers angle resolution", default=cs.GRIDRES)
    parser.add_argument("--precision", action='store_true', help="Use a more agressive control file to achieve better convergence")
    parser.add_argument("--test", action='store_true', help="Run a fast PELE_Platform test")
    parser.add_argument("--user_center", "-c", nargs='+', type=float, help='center of the box', default=None)
    parser.add_argument("--user_radius", "-r", type=float,  help="Radius of the box", default=None)
    parser.add_argument("--folder", "-wf", type=str,  help="Folder to apply the restart to", default=None)
    parser.add_argument("--pdb", action='store_true',  help="Use pdb files as output")
    parser.add_argument("--hbond", nargs='+',  help="Definition of kinase hbond", default= [None, None] )
    parser.add_argument("--msm", action="store_true",  help="Launch MSM")
    parser.add_argument("--out_in", action="store_true",  help="Launch outside inside adaptive")
    parser.add_argument("--induce_fit", action="store_true",  help="Launch induce fit adaptive")
    parser.add_argument("--adaptive", type=str,  help="Adaptive control_file")
    parser.add_argument("--pele", type=str,  help="Pele control_file")
    parser.add_argument("--precision_glide", type=str,  help="Glide precision.. Options = [SP, XP]", default="SP")
    parser.add_argument("--template", type=str,  help="External template for ligand", default=None)
    parser.add_argument("--rotamers", type=str,  help="External romtamers library for ligand", default=None)
    parser.add_argument("--lagtime", type=int,  help="MSM Lagtime", default=100)
    parser.add_argument("--msm_clust", type=int,  help="MSM cluster number", default=200)
    return parser.parse_args()


if __name__ == "__main__":
    args = parseargs()
    if args.hbond[0]:
        setattr(args, "software", "glide")
    elif args.adaptive and args.pele:
         setattr(args, "software", "adaptive")
    elif args.out_in:
        setattr(args, "software", "out_in")
    elif args.induce_fit:
        setattr(args, "software", "induce_fit")
    elif args.msm:
        setattr(args, "software", "msm")
    else:
        raise("Not specified action. Choose an option between msm/adaptive/out_in/induce_fit/hbond")
    platform_object = Launcher(args)
    platform_object.launch()
