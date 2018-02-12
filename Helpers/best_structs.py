import os
import argparse
import pandas as pd
import glob

"""

   Description: Parse all the reports found under 'path' and sort them all
   by the chosen criteria (Binding Energy as default) having into account the
   frequency our pele control file writes a structure through the -ofreq param
   (1 by default). To sort from higher to lower value use -f "max" otherwise
   will rank the structures from lower to higher criteria's values. The number
   of structures will be ranked is controlled by -i 'nstruct' (default 10).

   For any problem do not hesitate to contact us through the email address written below.

"""

__author__ = "Daniel Soler Viladrich"
__email__ = "daniel.soler@nostrumbiodiscovery.com"

# DEFAULT VALUES
ORDER = "min"
CRITERIA = "Binding Energy"
OUTPUT = "Structure_{}.pdb"
N_STRUCTS = 10
FREQ = 1
REPORT = "report"
TRAJ = "trajectory"
ACCEPTED_STEPS = 'numberOfAcceptedPeleSteps'
PATH = 'path'


def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument("path", type=str, help="Path to Pele's results root folder (Adaptive: path=/Pele/results/ Pele: path=/Pele/)")
    parser.add_argument("--crit", "-c", type=str, help="Criteria we want to rank and output the strutures for", default= CRITERIA)
    parser.add_argument("--nst", "-n", type=int, help="Number of produced structures" , default=N_STRUCTS)
    parser.add_argument("--sort", "-s", type=str, help="Look for minimum or maximum value --> Options: [min/max]", default=ORDER)
    parser.add_argument("--ofreq", "-f", type=int, help="Every how many steps the trajectory were outputted on PELE", default=FREQ)
    args = parser.parse_args()

    return args.path, args.crit, args.nst, args.sort, args.ofreq


def main(path, criteria="sasaLig", n_structs=500, sort_order="max", out_freq=FREQ):
    """

      Description: Rank the traj found in the report files under path
      by the chosen criteria. Finally, output the best n_structs.

      Input:

         Path: Path to look for *report* files in all its subfolders.

         Criteria: Criteria to sort the structures.
         Needs to be the name of one of the Pele's report file column.
         (Default= "Binding Energy")

         n_structs: Numbers of structures to create.

         sort_order: "min" if sorting from lower to higher "max" from high to low.

         out_freq: "Output frequency of our Pele control file"

     Output:

        f_out: Name of the n outpu
    """

    # Get Files
    all_reports = glob.glob(os.path.join(path, "*/*report*"))
    reports = [report for report in all_reports if(os.path.basename(os.path.dirname(report)).isdigit())]

    # Data Mining
    min_values = parse_values(reports, n_structs, criteria, sort_order)
    values = min_values[criteria].tolist()
    paths = min_values[PATH].tolist()
    epochs = [os.path.basename(os.path.normpath(os.path.dirname(Path))) for Path in paths]
    reports_indexes = min_values.report.tolist()
    step_indexes = min_values[ACCEPTED_STEPS].tolist()
    max_sasa_info = {i: [epoch, report, value, int(step)] for i, (epoch, report, value, step) in enumerate(zip(epochs, reports_indexes, values, step_indexes))}
    print(max_sasa_info)
    return max_sasa_info


def parse_values(reports, n_structs, criteria, sort_order):
    """

       Description: Parse the 'reports' and create a sorted array
       of size n_structs following the criteria chosen by the user.

    """

    INITIAL_DATA = [(PATH, []),
                    (REPORT, []),
                    (ACCEPTED_STEPS, []),
                    (criteria, [])
                    ]

    values = pd.DataFrame.from_items(INITIAL_DATA)
    for file in reports:
        report_number = os.path.basename(file).split("_")[-1]
        data = pd.read_csv(file, sep='    ', engine='python')
        selected_data = data.loc[:, [ACCEPTED_STEPS, criteria]]
        report_values = selected_data.nlargest(n_structs, criteria)
        report_values.insert(0, PATH, [file]*report_values[criteria].size)
        report_values.insert(1, REPORT, [report_number]*report_values[criteria].size)
        report_values = report_values[report_values[criteria].between(0.9, 1, inclusive=True)]
        try:
            values = pd.concat([values, report_values])
        except ValueError:
            values = report_values
    values.sort_values(criteria, ascending=False)
    return values


if __name__ == "__main__":
    path, criteria, interval, sort_order, out_freq = parse_args()
    main(path, criteria, interval, sort_order, out_freq)
