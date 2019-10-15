from argparse import ArgumentParser
import re
import prody
import os
import glob

def split_pdb(pdb_filename, output_path, prefix):
    struct = prody.parsePDB(pdb_filename)
    if prefix:
        out_pattern = "{}{}_{}"
    else:
        out_pattern = "{}{}{}"
    for i in range(struct.numCoordsets()):
        pattern = re.search(r"\w+/(\d+)/(\w+)_trajectory_(\d+)\.pdb", pdb_filename)
        # pattern = re.search(r"\w+/(\d+)/trajectory_(\d+)\.pdb", pdb_filename)
        if pattern is None:
            print("The file {} doesn't match the naming convention".format(pdb_filename))
            print("It won't be processed.")
            continue
        else:
            system_id = "{0[1]}_epoch_{0[0]}_processor_{0[2]}_model_{1}".format(pattern.groups(), i+1)
            # system_id = "epoch_{0[0]}_processor_{0[1]}_model_{1}".format(pattern.groups(), i + 1)
        outputfilename = out_pattern.format(output_path, prefix, system_id)
        prody.writePDB(outputfilename, struct, csets=i)


def main(input_folder, output_path, output_prefix):
    try: 
        os.path.isdir(output_path)
    except IOError:
        os.mkdir(output_path)
    for f in input_folder:
        if not os.path.isdir(f):
            print("The input {} isn't a folder it won't be used.".format(f))
            continue
        files_to_read = glob.glob(f + "*/*trajectory*.pdb")
        if len(files_to_read) == 0:
            files_to_read = glob.glob(f + "*/*trajectory*.xtc")
            if len(files_to_read) == 0:
                print("No file in a valid format (.pdn or .xtc) found. \nDiscontinuing the folder {}".format(f))
            else:
                print("No treatment for .xtc files implemented yet")
        else:
            for pdb_file in files_to_read:
                split_pdb(pdb_file, output_path, output_prefix)

if __name__ == '__main__':
    parser = ArgumentParser(description="A small script to split ALL the trajectories files into files.")
    parser.add_argument("input_folder", help="The folder to process", nargs='+')
    parser.add_argument("-output_folder", default="", help="The folder where the sturctures should be written.")
    parser.add_argument("-output_prefix", default="", help="A prefix for the name of the written sturctures.")
    args = parser.parse_args()
    main(args.input_folder, args.output_folder, args.output_prefix)
