import matplotlib
matplotlib.use("Agg")
import argparse
import os
import pele_platform.MSM.main as msm
import pele_platform.constants as cs
import pele_platform.Rescore.main as gl
import pele_platform.Rescore.simulation as ad


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

        elif args.software in ["glide", "induce_fit"]:
            gl.run(args)

        elif args.software == "frag":
            main = os.path.join(cs.DIR, "LigandGrowing/grow.py")
            "{} {} -cp {} -fp {} -ca {} -fa {}".format(cs.PYTHON3, main, args.system, args.frag, args.ca, args.fa)



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
    parser.add_argument('--confile', type=str, help='your own pele configuration file', default=None)
    parser.add_argument('--adapt_conf', type=str, help='your own adaptive pele configuration file', default=None)
    parser.add_argument('--native', type=str, help='native file to compare RMSD to', default="")
    parser.add_argument('--cpus', type=int, help='number of processors', default=cs.CPUS)
    parser.add_argument("--core", type=int, help="Give one atom of the core section", default=-1)
    parser.add_argument("--mtor", type=int, help="Gives the maximum number of torsions allowed in each group.  Will freeze bonds to extend the core if necessary.", default=4)
    parser.add_argument("--n", type=int, help="Maximum Number of Entries in Rotamer File", default=1000)
    parser.add_argument("--clean", help="Whether to clean up all the intermediate files", action='store_true')
    parser.add_argument("--restart", type=str, help="Restart the platform from [all, pele, msm] with these keywords", default=cs.PLATFORM_RESTART)
    parser.add_argument("--gridres", type=str, help="Rotamers angle resolution", default=cs.GRIDRES)
    parser.add_argument("--precision", action='store_true', help="Use a more agressive control file to achieve better convergence")
    parser.add_argument("--test", action='store_true', help="Run a fast pele_platform test")
    parser.add_argument("--user_center", "-c", nargs='+', type=str, help='center of the box', default=None)
    parser.add_argument("--box_radius", "-r", type=float,  help="Radius of the box", default=None)
    parser.add_argument("--folder", "-wf", type=str,  help="Folder to apply the restart to", default=None)
    parser.add_argument("--pdb", action='store_true',  help="Use pdb files as output")
    parser.add_argument("--hbond", nargs='+',  help="Definition of kinase hbond", default= [None, None] )
    parser.add_argument("--msm", action="store_true",  help="Launch MSM")
    parser.add_argument("--out_in", action="store_true",  help="Launch outside inside adaptive")
    parser.add_argument("--full", action="store_true",  help="Launch ful ful binding site exploration")
    parser.add_argument("--in_out", action="store_true",  help="Launch inside outside soft adaptive")
    parser.add_argument("--in_out_soft", action="store_true",  help="Launch inside outside adaptive")
    parser.add_argument("--water_exp", type=str,  help="Launch water exploration adaptive PELE", default=None)
    parser.add_argument("--water_lig", nargs="+",  help="Launch ligand-water exploration adaptive PELE", default=None)
    parser.add_argument("--water_center", nargs="+",  help="Launch ligand-water exploration adaptive PELE", default=None)
    parser.add_argument("--induce_fit", action="store_true",  help="Launch induce fit adaptive")
    parser.add_argument("--adaptive", type=str,  help="Adaptive control_file")
    parser.add_argument("--pele", type=str,  help="Pele control_file")
    parser.add_argument("--precision_glide", type=str,  help="Glide precision.. Options = [SP, XP]", default="SP")
    parser.add_argument("--template", type=str,  help="External template for ligand", default=None)
    parser.add_argument("--rotamers", type=str,  help="External romtamers library for ligand", default=None)
    parser.add_argument("--lagtime", type=int,  help="MSM Lagtime", default=100)
    parser.add_argument("--msm_clust", type=int,  help="MSM cluster number", default=200)
    parser.add_argument("--frag", type=str,  help="Fragment pdb")
    parser.add_argument("--ca", type=str,  help="Core Atom")
    parser.add_argument("--fa", type=str,  help="Fragment Atom")
    parser.add_argument("--noeq", action="store_false",  help="Whether to do a initial equilibration or not")
    parser.add_argument("--skip_prep", action="store_true",  help="Whether to do the initial preprocessing or not")
    parser.add_argument("--nonstandard", nargs="+",  help="Mid Chain non standard residues to be treated as ATOM not HETATOM", default = [])
    parser.add_argument('--solvent', type=str, help='Type of implicit solvent (OBC/VDGBNP). default [OBC]. i.e. --solvent VDGBNP', default="OBC")
    parser.add_argument('--randomize', action="store_true", help='Randomize ligand position around protein')
    parser.add_argument("--atom_dist", nargs="+",  help="Number of the atoms to calculate the distance in between i.e --atom dist 123 456", default=None)
    parser.add_argument('--input', nargs="+", help='Set initial input for simulation')

    return parser.parse_args()

def set_software_to_use():
    """
    Auxiliar Function to set low variable software
    which will be use to handle differences 
    betwwen PELE features along the program
    """
    if args.hbond[0]:
        setattr(args, "software", "glide")
    elif args.out_in or args.water_lig or args.full or args.water_exp or args.in_out_soft or args.in_out or args.induce_fit or  (args.adaptive and args.pele):
        setattr(args, "software", "adaptive")
    elif args.msm:
        setattr(args, "software", "msm")
    elif args.frag and args.ca and args.fa:
        setattr(args, "software", "frag")
    else:
        raise ValueError("Not specified action. Choose an option between msm/adaptive/out_in/induce_fit/hbond")


if __name__ == "__main__":
    args = parseargs()
    set_software_to_use()
    Launcher(args).launch()
