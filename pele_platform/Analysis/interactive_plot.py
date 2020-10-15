import argparse
import AdaptivePELE.analysis.interactivePlot as plot


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--x", type=int, required=True, help="x-axis values", default=5)
    parser.add_argument("--y", type=int, required=True, help="y-axis values", default=6)
    parser.add_argument("--steps", type=int, required=False, help="Adaptive steps", default=8)
    parser.add_argument("--output", type=str, required=False, help="Output folder")


if __name__ == "__main__":
    x, y, z, output = parse_args()
    plot.main(x, y, z, path=current_dir, out_freq=1, output=output_folder, numfolders=False, topology=None)
