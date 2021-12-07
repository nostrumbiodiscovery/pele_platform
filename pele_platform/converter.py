# -*- coding: utf-8 -*-
"""
This is the converter module of the platform that can take PELE trajectories
and convert them to different trajectory formats.
"""

__author__ = "Nostrum Biodiscovery"
__email__ = "pelesupport@nostrumbiodiscovery.com"
__license__ = "Apache-2.0"
__version__ = "1.0"


# Constants
FORMATS = ('pdb', 'xtc')


# Methods
def parse_args():
    """
    Command line parser.

    Returns
    -------
    input_path : str
        Path with PELE trajectories to convert
    output_path : str
        Path where converted trajectories will be saved
    input_format : str
        Original format to convert
    output_format : str
        New format that is requested
    topology_path : str
        Path to topology PDB file
    delete : bool
        Whether to delete original trajectory after conversion or not
    verify : bool
        Whether to verify new trajectories after conversion or not
    n_processors : int
        Number of processors to use during the conversion
    trajectory_name : str
        Name of PELE's trajectory files
    """

    # Parser setup
    from argparse import ArgumentParser

    parser = ArgumentParser(
        description='PELE Platform Converter')

    parser.add_argument('input_path', type=str,
                        help='Path with PELE trajectories to convert')
    parser.add_argument('-o', '--output_path', type=str, default=None,
                        help='Path where converted trajectories will be saved')
    parser.add_argument('-if', '--input_format', type=str, required=True,
                        choices=FORMATS,
                        help='Original format to convert')
    parser.add_argument('-of', '--output_format', type=str, required=True,
                        choices=FORMATS,
                        help='New format that is requested')
    parser.add_argument('-t', '--topology', type=str, default=None,
                        help='Path to topology PDB file')
    parser.add_argument('-d', '--delete', dest='delete',
                        action='store_true',
                        help='Delete original trajectory after conversion')
    parser.add_argument('-n', '--n_processors', type=int,
                        default=1,
                        help='Number of processors to use during the ' +
                        'conversion')
    parser.add_argument('--trajectory_name', type=str,
                        default='trajectory',
                        help='Name of PELE\'s trajectory files')

    group = parser.add_mutually_exclusive_group()
    group.add_argument('--verify', dest='verify',
                        action='store_true',
                        help='New trajectories are double checked for ' +
                             'consistency after conversion')
    group.add_argument('--dont_verify', dest='verify',
                        action='store_false',
                        help='New trajectories are double checked for ' +
                             'consistency after conversion')

    # Set defaults
    parser.set_defaults(delete=False, verify=None)

    # Parse arguments
    parsed_args = parser.parse_args()

    # Extract arguments
    input_path = parsed_args.input_path
    output_path = parsed_args.output_path
    input_format = parsed_args.input_format
    output_format = parsed_args.output_format
    topology_path = parsed_args.topology
    delete = parsed_args.delete
    n_processors = parsed_args.n_processors
    trajectory_name = parsed_args.trajectory_name

    # Default is to only verify when deleting original trajectories
    if parsed_args.verify is None:
        verify = delete
    else:
        verify = parsed_args.verify

    return input_path, output_path, input_format, output_format, \
        topology_path, delete, verify, n_processors, trajectory_name


def print_args(input_path, output_path, input_format, output_format,
               topology_path, delete, verify, n_processors,
               trajectory_name):
    """
    It prints the parameters supplied by the user.

    Parameters
    ----------
    input_path : str
        Path with PELE trajectories to convert
    output_path : str
        Path where converted trajectories will be saved
    input_format : str
        Original format to convert
    output_format : str
        New format that is requested
    topology_path : str
        Path to topology PDB file
    delete : bool
        Whether to delete original trajectory after conversion or not
    verify : bool
        Whether to verify new trajectories after conversion or not
    n_processors : int
        Number of processors to use during the conversion
    trajectory_name : str
        Name of PELE's trajectory files
    """
    max_len = max([len(str(item)) if item is not None else 1
                   for item in [input_path, output_path,
                                input_format, output_format,
                                topology_path, delete,
                                verify, n_processors,
                                trajectory_name]])

    print('-> Input parameters:')
    print(' - input_path:      ',
          ' ' * (max_len - (len(str(input_path))
                            if input_path is not None else 1)),
          '-' if input_path is None else input_path)
    print(' - output_path:     ',
          ' ' * (max_len - (len(str(output_path))
                            if output_path is not None else 1)),
          '-' if output_path is None else output_path)
    print(' - input_format:    ',
          ' ' * (max_len - (len(str(input_format))
                            if input_format is not None else 1)),
          '-' if input_format is None else input_format)
    print(' - output_format:   ',
          ' ' * (max_len - (len(str(output_format))
                            if output_format is not None else 1)),
          '-' if output_format is None else output_format)
    print(' - topology_path:   ',
          ' ' * (max_len - (len(str(topology_path))
                            if topology_path is not None else 1)),
          '-' if topology_path is None else topology_path)
    print(' - delete:          ',
          ' ' * (max_len - (len(str(delete))
                            if delete is not None else 1)),
          '-' if delete is None else delete)
    print(' - verify:          ',
          ' ' * (max_len - (len(str(verify))
                            if verify is not None else 1)),
          '-' if verify is None else verify)
    print(' - n_processors:    ',
          ' ' * (max_len - (len(str(n_processors))
                            if n_processors is not None else 1)),
          '-' if n_processors is None else n_processors)
    print(' - trajectory_name: ',
          ' ' * (max_len - (len(str(trajectory_name))
                            if trajectory_name is not None else 1)),
          '-' if trajectory_name is None else trajectory_name)
    print()


def convert_trajectories(input_path, output_path, input_format,
                         output_format, topology_path, n_processors,
                         trajectory_name, verify):
    """
    It converts trajectories to the chosen format.

    Parameters
    ----------
    input_path : str
        Path with PELE trajectories to convert
    output_path : str
        Path where converted trajectories will be saved
    input_format : str
        Original format to convert
    output_format : str
        New format that is requested
    topology_path : str
        Path to topology PDB file
    n_processors : int
        Number of processors to use during the conversion
    trajectory_name : str
        Name of PELE's trajectory files
    verify : bool
        Verify that written trajectory is correctly built and matches with
        original one

    Returns
    -------
    input_sizes : int
        Total size of input trajectories
    output_sizes : int
        Total size of output trajectories
    converted_trajectories : list[str]
        List of original trajectories that have been converted successfully
    """
    import os
    import glob
    from functools import partial
    from multiprocessing import Pool
    from tqdm import tqdm

    # Initialize primary variables
    epoch_dirs = glob.glob(os.path.join(input_path, '[0-9]*'))

    # Filter out non digit folders
    epochs = [os.path.basename(path) for path in epoch_dirs
              if os.path.basename(path).isdigit() and os.path.isdir(path)]

    # Sort epochs by number
    epochs = sorted(epochs, key=int)

    # Tweak to read a directory from standard PELE (not coming
    # from adaptive)
    if len(epochs) == 0:
        epochs = ['']

    if output_path is None:
        general_output_path = input_path
        if os.path.isfile(general_output_path):
            general_output_path = os.path.dirname(input_path)
    else:
        general_output_path = output_path

    input_sizes = 0
    output_sizes = 0
    converted_trajectories = list()
    messages = list()  # store all warning messages and print them at the end

    for epoch in tqdm(epochs, desc=" - Progress: "):
        if epoch == '':
            trajectory_path = input_path
        else:
            trajectory_path = os.path.join(input_path, epoch)

        # Case of a single trajectory in a single directory
        if os.path.isfile(trajectory_path):
            trajectory_prefix = \
                os.path.splitext(os.path.basename(input_path))[0]
            specific_output_path = os.path.join(general_output_path,
                                                trajectory_prefix + '.' +
                                                output_format)
            input_size, output_size = \
                convert_trajectory(trajectory_path, specific_output_path,
                                   topology_file=topology_path,
                                   verify=verify)

            input_sizes += input_size
            output_sizes += output_size
            if output_size > 0:
                converted_trajectories.append(trajectory_path)
            else:
                messages.append(f'Warning, verification failed for ' +
                                f'{trajectory_path}. It could not ' +
                                f'be converted properly and original ' +
                                f'trajectory will not be deleted')
            break

        # Case of a single trajectory in multiple Adaptive directories
        trajectory_path = os.path.join(trajectory_path,
                                       trajectory_name + '.' + input_format)
        specific_output_path = os.path.join(general_output_path, epoch,
                                            trajectory_name + '.' +
                                            output_format)

        if os.path.isfile(trajectory_path):
            input_size, output_size = \
                convert_trajectory(trajectory_path, specific_output_path,
                                   topology_file=topology_path,
                                   verify=verify)

            input_sizes += input_size
            output_sizes += output_size
            if output_size > 0:
                converted_trajectories.append(trajectory_path)
            else:
                messages.append(f'Warning, verification failed for ' +
                                f'{trajectory_path}. It could not ' +
                                f'be converted properly and original ' +
                                f'trajectory will not be deleted')
            continue

        # Case of multiple trajectories in multiple Adaptive directories
        trajectory_path = os.path.join(input_path, epoch, trajectory_name)
        trajectories = glob.glob(trajectory_path + '_[0-9]*.' + input_format)

        if len(trajectories) == 0:
            messages.append(f'Warning, no {input_format.upper()} ' +
                            f'trajectories were found ' +
                            f'at {os.path.dirname(trajectory_path)}')

        parallel_function = partial(convert_trajectory_in_parallel,
                                    general_output_path, epoch,
                                    output_format, topology_path,
                                    verify)

        with Pool(n_processors) as pool:
            sizes = pool.map(parallel_function, trajectories)

        for (input_size, output_size), trajectory in zip(sizes,
                                                         trajectories):
            input_sizes += input_size
            output_sizes += output_size
            if output_size > 0:
                converted_trajectories.append(trajectory)
            else:
                messages.append(f'Warning, verification failed for ' +
                                f'{trajectory}. It could not be converted ' +
                                f'properly and original trajectory will ' +
                                f'not be deleted')

    # Print all warnings at the end
    for message in messages:
        print(message)

    return input_sizes, output_sizes, converted_trajectories


def convert_trajectory_in_parallel(general_output_path, epoch,
                                   output_format, topology_path,
                                   verify, trajectory_path):
    """
    It converts trajectories to the chosen format, in parallel.

    Parameters
    ----------
    general_output_path : str
        General path where converted trajectories will be saved
    epoch : int
        Epoch the current trajectory belongs to
    output_format : str
        New format that is requested
    topology_path : str
        Path to topology PDB file
    verify : bool
        Verify that written trajectory is correctly built and matches with
        original one
    trajectory_path : str
        Path of PELE's trajectory to convert

    Returns
    -------
    input_size : int
        Size of input trajectory
    output_size : int
        Size of output trajectory
    """
    import os

    trajectory_prefix = os.path.splitext(os.path.basename(trajectory_path))[0]
    specific_output_path = os.path.join(general_output_path, str(epoch),
                                        trajectory_prefix + '.' +
                                        output_format)

    if os.path.isfile(trajectory_path):
        input_size, output_size = \
            convert_trajectory(trajectory_path, specific_output_path,
                               topology_file=topology_path,
                               verify=verify)

    return input_size, output_size


def convert_trajectory(trajectory_path, output_path, topology_file=None,
                       verify=True):
    """
    It converts the trajectory to the format defined in the output path.

    Parameters
    ----------
    trajectory_path : str
        Path to the trajectory to convert
    output_path : str
        Output path of the trajectory, including the format extension
        that is requested
    topology_file : str
        Path to topology file. Default is None
    verify : bool
        Verify that written trajectory is correctly built and matches with
        original one

    Returns
    -------
    input_size : int
        Size of input trajectory
    output_size : int
        Size of output trajectory
    """
    import os
    import mdtraj

    # Check path
    parent_path = os.path.dirname(output_path)

    if not os.path.isdir(parent_path):
        # Workaround to create missing folders in parallel
        try:
            os.makedirs(parent_path)
        except FileExistsError:
            pass

    # Check extension
    extension1 = \
        os.path.splitext(os.path.basename(trajectory_path))[-1].strip('.')
    if extension1 not in FORMATS:
        raise ValueError('Invalid input trajectory format. ' +
                         f'Expected one of {str(FORMATS)}, got {extension1}')
    extension2 = \
        os.path.splitext(os.path.basename(output_path))[-1].strip('.')
    if extension2 not in FORMATS:
        raise ValueError('Invalid output trajectory format. ' +
                         f'Expected one of {str(FORMATS)}, got {extension2}')

    if topology_file is not None and extension1 == 'xtc':
        trajectory = mdtraj.load(trajectory_path, top=topology_file)
    else:
        trajectory = mdtraj.load(trajectory_path)

    trajectory.save(output_path)

    # Get sizes
    input_size = os.path.getsize(trajectory_path)
    output_size = os.path.getsize(output_path)

    # Verify new trajectory
    if verify:
        is_okay = verify_trajectory(trajectory_path, output_path)
    else:
        is_okay = True

    if is_okay:
        return input_size, output_size
    else:
        if os.path.isfile(output_path):
            os.remove(output_path)
        return input_size, 0


def verify_trajectory(trajectory_path, output_path):
    """
    It verifies that written trajectory contains the same information as
    the original one.

    Parameters
    ----------
    trajectory_path : str
        Path to the original trajectory
    output_path : str
        Path to the new trajectory

    Returns
    -------
    is_okay : bool
        Whether new built trajectory is correct or not
    """
    import os
    import mdtraj
    import numpy as np

    # One of them must be a PDB file.
    extension1 = os.path.splitext(os.path.basename(
        trajectory_path))[-1].strip('.')
    extension2 = os.path.splitext(os.path.basename(
        output_path))[-1].strip('.')

    if extension1 == 'pdb':
        pdb_trajectory = trajectory_path
        xtc_trajectory = output_path
    elif extension2 == 'pdb':
        pdb_trajectory = output_path
        xtc_trajectory = trajectory_path
    else:
        raise ValueError('Either input format or output format must be PDB')

    trajectory1 = mdtraj.load(pdb_trajectory)
    trajectory2 = mdtraj.load(xtc_trajectory, top=trajectory1)

    is_okay = np.allclose(trajectory1.xyz, trajectory2.xyz, atol=1.e-3)

    return is_okay


def delete_trajectories(trajectories_to_delete):
    """
    It deletes original trajectories.

    Parameters
    ----------
    trajectories_to_delete : list[str]
        List of trajectories to delete

    Returns
    -------
    deleted_size : int
        Total size of deleted trajectories
    """
    import os

    deleted_size = 0

    for trajectory in trajectories_to_delete:
        if os.path.isfile(trajectory):
            deleted_size += os.path.getsize(trajectory)
            os.remove(trajectory)
        else:
            raise ValueError(f'Trajectory {trajectory} was not found and ' +
                             f'could not be deleted')

    return deleted_size


def bytes_to_string(bytes):
    """
    It generates a string with a proper format to represent bytes.

    Parameters
    ----------
    bytes : int
        A quantity of bytes

    Returns
    -------
    size_str : str
        The string representing the number of bytes with a proper format
    """
    kilobytes = bytes / 1000

    if kilobytes < 1.0:
        return '{:.2f} B'.format(bytes)

    megabytes = kilobytes / 1000.0
    if megabytes < 1.0:
        return '{:.2f} kB'.format(kilobytes)

    gigabytes = megabytes / 1000
    if gigabytes < 1.0:
        return '{:.2f} MB'.format(megabytes)

    terabytes = gigabytes / 1000
    if terabytes < 1.0:
        return '{:.2f} GB'.format(gigabytes)

    return '{:.2f} TB'.format(terabytes)


# Main workflow to be executed
if __name__ == "__main__":
    # Import external libraries
    import os

    # Print header
    from pele_platform.constants import constants
    print(constants.converter_version_header.format(__version__))

    input_path, output_path, input_format, output_format, topology_path, \
        delete, verify, n_processors, trajectory_name = parse_args()

    print_args(input_path, output_path, input_format, output_format,
               topology_path, delete, verify, n_processors,
               trajectory_name)

    # Check input and output paths
    if not os.path.exists(input_path):
        raise ValueError('Invalid input path, it does not exist')

    if output_path is not None and not os.path.exists(output_path):
        os.makedirs(output_path)

    # Check input and output formats
    if input_format == output_format:
        raise ValueError('Input and output trajectory formats are equal. ' +
                         'Nothing to do.')

    # Check topology
    if (input_format == 'xtc' and
            (topology_path is None or not os.path.isfile(topology_path))):
        raise ValueError('XTC trajectories cannot be read without a ' +
                         'topology file. It must be supplied.')

    # Convert trajectories
    print('-> Converting trajectories:')
    input_size, output_size, converted_trajectories = \
        convert_trajectories(input_path, output_path, input_format,
                             output_format, topology_path,
                             n_processors, trajectory_name, verify)

    print(' - Read ' + bytes_to_string(input_size))
    print(' - Written ' + bytes_to_string(output_size))
    print()

    # Delete original trajectories
    if delete:
        print('-> Deleting original trajectories:')
        deleted_size = delete_trajectories(converted_trajectories)
        print(f' - Deleted ' + bytes_to_string(deleted_size))
        print()
