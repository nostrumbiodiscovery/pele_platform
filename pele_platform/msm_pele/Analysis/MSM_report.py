from fpdf import FPDF
import re
import sys
import os
import subprocess
import msm_pele.Helpers.helpers as hp
import argparse
from AdaptivePELE.analysis import plotAdaptive

CONTACTS = "contacts"
steps_Run, Xcol, Ycol, filename, kind_Print, colModifier, traj_range = 8, 6, 5, "report_", "PRINT_RMSD_STEPS", 4, None


def arg_parse():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('control_file', type=str, help='adaptive control file')
    args = parser.parse_args()
    return args


class Plots():
    
    def __init__(self, folder):
        self.pmf = os.path.join(folder, "PMF_plots/pmf_run_0.png")
        self.probability = os.path.join(folder, "eigenvector_plots/eigenvector_1_alone_run_0.png")
        self.transitions = os.path.join(folder, "transitions/transitions.png")
        self.transition_hist = os.path.join(folder, "transitions/transition_hist.png")
        


def report_MSM(env, folder):


    OUTPUT = os.path.join(env.results,  'MSM_report.pdf')
    dG = get_dG_line(env)
    plots = Plots(folder)

    pdf = FPDF()
    pdf.add_page() 

    # Title
    pdf.set_font('Arial', 'B', 15)
    pdf.cell(80)
    pdf.cell(30, 10, "{} {} kcal/mol".format(env.residue, dG), align='C')

    # PMF
    pdf.cell(-74)
    pdf.set_font('Arial', 'B', 11)
    pdf.cell(10, 49, 'PMF' +70*"\t" + "Probability")
    pdf.image(plots.pmf, 10, 40, 83)

    # Probability
    pdf.cell(0)
    pdf.set_font('Arial', 'B', 11)
    pdf.image(plots.probability, 100, 40, 83)

    # Transitions
    pdf.cell(-170)
    pdf.set_font('Arial', 'B', 11)
    pdf.cell(10, 249, 'Transition Hist')
    pdf.cell(0)
    pdf.set_font('Arial', 'B', 11)
    pdf.image(plots.transition_hist, 10, 140, 83)

    #Output report    
    pdf.output(OUTPUT, 'F')


def get_dG_line(env):
    with open(os.path.join(env.results, "results.txt"), "r") as f:
        lines = f.readlines()
        for line in reversed(lines):
            if not line.startswith("#"):
                try:
                    _, dg, std, _, _, _ = line.strip("\n").split()
                except ValueError:
                    _, dg, std, _, _ = line.strip("\n").split()
                return dg + " +- " + std 

    
def plot(Xcol, Ycol, path, name, xcol_name, ycol_name, zcol=5, style="PRINT_BE_RMSD"):
    plot_line = plotAdaptive.generatePrintString(1000, Xcol, Ycol, "report_", style, zcol, None).strip("\n")
    if "$" in plot_line:
        plot_line = plot_line.replace("$", "$ ")
    command = '''gnuplot -e "set terminal png; set output '{}'; set xlabel '{}'; set ylabel '{}'; {}"'''.format(
        name, xcol_name, ycol_name, plot_line)
    os.system(command)


if __name__ == "__main__":
    args = arg_parse()
    main(args.control_file)
