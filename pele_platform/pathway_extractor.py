# -*- coding: utf-8 -*-
"""
This is the pathway extractor module of the platform that can extract
ligand pathways from a PELE simulation.
"""

__author__ = "Nostrum Biodiscovery"
__email__ = "pelesupport@nostrumbiodiscovery.com"
__license__ = "Apache-2.0"
__version__ = "1.0"


# Methods
def parse_args():
    """
    Command line parser.

    Returns
    -------
    epoch : str
        Path to the epoch to search the snapshot
    trajectory : int
        Trajectory number of the snapshot to extract
    snapshot : int
        Snapshot to select (in accepted steps)
    output : str
        Output path where to write the resulting pathway
    name : str
        Name of the PDB to write the resulting pathway
    topology : str
        Path to topology PDB file for loading non-PDB trajectories
    """

    # Parser setup
    from argparse import ArgumentParser

    parser = ArgumentParser(
        description='PELE Platform Pathway Extractor')

    parser.add_argument("epoch", type=str,
                        help="Path to the epoch to search the snapshot")
    parser.add_argument("trajectory", type=int,
                        help="Trajectory number of the snapshot to extract")
    parser.add_argument("snapshot", type=int,
                        help="Snapshot to select (in accepted steps)")
    parser.add_argument("-o", type=str, default=None,
                        help="Output path where to write the resulting " +
                             "pathway")
    parser.add_argument("--name", type=str, default="pathway.pdb",
                        help="Name of the PDB to write the resulting" +
                             "pathway")
    parser.add_argument("--top", type=str, default=None,
                        help="Name of the PDB topology for loading " +
                             "non-PDB trajectories")

    args = parser.parse_args()

    trajectory = args.trajectory
    snapshot = args.snapshot
    epoch = args.epoch
    output = args.o
    name = args.name
    topology = args.top

    return trajectory, snapshot, epoch, output, name, topology


def print_args(trajectory, snapshot, epoch, output, name, topology):
    """
    It prints the parameters supplied by the user.

    Parameters
    ----------
    epoch : str
        Path to the epoch to search the snapshot
    trajectory : int
        Trajectory number of the snapshot to extract
    snapshot : int
        Snapshot to select (in accepted steps)
    output : str
        Output path where to write the resulting pathway
    name : str
        Name of the PDB to write the resulting pathway
    topology : str
        Path to topology PDB file for loading non-PDB trajectories
    """
    max_len = max([len(str(item)) if item is not None else 1
                   for item in [epoch, trajectory,
                                snapshot, output,
                                name, topology]])

    print('-> Input parameters:')
    print(' - epoch:             ',
          ' ' * (max_len - (len(str(epoch))
                            if epoch is not None else 1)),
          '-' if epoch is None else epoch)
    print(' - trajectory:        ',
          ' ' * (max_len - (len(str(trajectory))
                            if trajectory is not None else 1)),
          '-' if trajectory is None else trajectory)
    print(' - snapshot:          ',
          ' ' * (max_len - (len(str(snapshot))
                            if snapshot is not None else 1)),
          '-' if snapshot is None else snapshot)
    print(' - output:            ',
          ' ' * (max_len - (len(str(output))
                            if output is not None else 1)),
          '-' if output is None else output)
    print(' - name:              ',
          ' ' * (max_len - (len(str(name))
                            if name is not None else 1)),
          '-' if name is None else name)
    print(' - topology:          ',
          ' ' * (max_len - (len(str(topology))
                            if topology is not None else 1)),
          '-' if topology is None else topology)
    print()


# Main workflow to be executed
if __name__ == "__main__":
    # Import external libraries
    import os

    # Print header
    from pele_platform.constants import constants
    print(constants.pathway_extractor_version_header.format(__version__))

    # Read CLI arguments
    trajectory, snapshot, epoch, output, name, topology = parse_args()

    # Print them
    print_args(trajectory, snapshot, epoch, output, name, topology)

    # Check input and output paths
    if not os.path.exists(epoch):
        raise ValueError('Invalid epoch path, it does not exist')

    if output is not None and not os.path.exists(output):
        os.makedirs(output)

    # Run Adaptive script
    from AdaptivePELE.analysis.backtrackAdaptiveTrajectory \
        import main as backtrack_trajectory

    try:
        backtrack_trajectory(trajectory, snapshot, epoch, output,
                             name, topology)
    except TypeError:
        raise TypeError('Trajectory generator failed. Remember to supply a '
                        'topology file when reading XTC trajectories')
