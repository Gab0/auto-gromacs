import argparse


def parse_arguments():
    """Parse commandline arguments."""
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-f",
        "--file-prefix",
        dest='FilePrefix',
        nargs="*",
        help="File prefix containing multiple GROMACS simulation directories."
    )

    parser.add_argument(
        "-d",
        "--auto-detect",
        dest='AutoDetect',
        nargs="*",
        help=""
    )

    parser.add_argument(
        "-t",
        "--traj-suffix",
        dest='TrajSuffix',
        default="",
        help=""
    )

    parser.add_argument(
        "-M",
        "--matrix-only",
        dest='matrix_only',
        action="store_true",
        help=""
    )

    parser.add_argument(
        "-T",
        "--timeseries",
        dest='DoTimeseries',
        action="store_true",
        help=""
    )

    parser.add_argument(
        "-w",
        "--write",
        dest='WriteOutput',
        action="store_true",
        help=""
    )

    parser.add_argument(
        "-i",
        "--identifier",
        dest='OutputIdentifier',
        required=True,
        help="Unique identifier for the current analysis." +
        " Will be included in all output filenames."
    )

    parser.add_argument(
        '-m',
        "--reference-mean",
        dest='ReferenceMean',
        action="store_true",
        help=""
    )

    parser.add_argument(
        '-s',
        "--selection",
        dest='SimulationSelection',
        help="Simulation directories to be included in the analysis." +
        " Will be selected interactively if not specified."

    )

    parser.add_argument(
        '--load-session',
        help="Filepath to load session files from (using pickle)."
    )

    parser.add_argument(
        '--write-session',
        help="Filepath to write session files to (using pickle)."
    )

    parser.add_argument(
        '--no-plot',
        action="store_true",
        help="When TRUE it will NOT create plots."
    )

    parser.add_argument(
        '--reference-structure',
        help="Path to reference .pdb structure to derive mutations from.\n" +
        "This allows the display of structure-specific mutations in output plots."
    )

    return parser.parse_args()
