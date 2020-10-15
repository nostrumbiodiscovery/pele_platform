import argparse
import os
import AdaptivePELE.analysis.interactivePlot as plot


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("x", type=int, help="x-axis values", default=5)
    parser.add_argument("y", type=int, help="y-axis values", default=6)
    parser.add_argument("--steps", type=int, required=False, help="Adaptive steps", default=8)
    parser.add_argument("--output", type=str, required=False, help="Output folder", default="interactive_output")
    parser.add_argument("--topology", type=str, required=False, help="Topology file", default=None)
    args = parser.parse_args()
    
    return args.x, args.y, args.steps, args.output, args.topology

if __name__ == "__main__":
    x, y, steps, output_folder, topology = parse_args()
    current_dir = os.getcwd()
    plot.main(x, y, steps, path=current_dir, out_freq=1, output=output_folder, numfolders=False, topology=topology)
