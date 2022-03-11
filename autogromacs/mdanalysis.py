from typing import List, Optional, cast, Tuple
import enum
import argparse
import sys
import os
import re
import warnings
import freesasa
import numpy as np

import MDAnalysis as mda
from MDAnalysis.analysis import align, rms, pca, psa

from . import mdplots

warnings.filterwarnings("ignore", category=DeprecationWarning)

STANDARD_SELECTION = "protein and name CA"


class AnalysisSession():
    """Holds what is saved from all loaded Universes for a single analysis."""
    labels: List[str] = []
    rmsd_series: List[np.ndarray] = []
    rmsf_series: List[np.ndarray] = []
    sample_rmsf: List[np.ndarray] = []
    pca_series: List[np.ndarray] = []
    total_times: List[int] = []
    samples: List[mda.Universe] = []
    sasa: List[np.ndarray] = []


class SeriesMode(enum.Enum):
    MONOLITHIC = 0
    STACKED = 1


class OperationMode():
    """
    Handles the operation mode for the anaylsis,
    which is based on the user input arguments.
    """
    compare_pairwise = True
    compare_timeseries = True

    def __init__(self, arguments):
        if arguments.matrix_only:
            self.compare_pairwise = True
            self.compare_timeseries = False



class AlignType(enum.Enum):
    """
    Extracts a single frame structural summary from
    trajectories by using different methods.
    """
    FIRST_FRAME = 0
    MEAN_FRAME = 1

    def extract(self, traj):
        """
        Extracts a single frame representation from a trajectory,
        based on the method represented by the instantiated Enum (self).
        """
        if self == AlignType.FIRST_FRAME:
            return traj[:, 0, :][:, None, :]
        if self == AlignType.MEAN_FRAME:
            return traj.mean(axis=1)[:, None, :]


def parse_arguments():
    """Parse commandline arguments."""
    parser = argparse.ArgumentParser()

    parser.add_argument("-f", dest='FilePrefix', nargs="*")
    parser.add_argument("-d", dest='AutoDetect', nargs="*")

    parser.add_argument("-t", dest='TrajSuffix', default="")
    parser.add_argument("-M", dest='matrix_only', action="store_true")
    parser.add_argument("-T", dest='DoTimeseries', action="store_true")

    parser.add_argument("-w", dest='WriteOutput', action="store_true")
    parser.add_argument("-i", "--identifier",
                        dest='OutputIdentifier', required=True)

    parser.add_argument('-m', dest='ReferenceMean', action="store_true")


    parser.add_argument('-s', dest='SimulationSelection')
    return parser.parse_args()


def autodetect_files(root_path, pattern="md.gro") -> List[str]:
    """Autodetect GROMACS simulation directories inside a 'project' folder."""
    def file_to_prefix(f):
        ext = f.split(".")[-1]
        return f.replace("." + ext, "")

    def detect(path):
        for f in os.listdir(path):
            F = os.path.join(path, f)
            if os.path.isdir(F):
                yield from detect(F)
            elif f.endswith(pattern):
                yield file_to_prefix(F)

    return list(sorted(detect(root_path)))


class Positions():
    def __init__(self, pos, universe_name, traj_idx):
        self.positions = pos
        self.universe_name = universe_name
        self.traj_idx = traj_idx

        self.label = self.make_label()

    def make_label(self):
        return "_".join([
            clean_universe_prefix(self.universe_name),
            index_label(self.traj_idx)
        ])


def index_label(label_id: int) -> str:
    common = {
        0: "initial",
        -1: "last"
    }
    try:
        return common[label_id]
    except KeyError:
        return str(label_id)


def clean_universe_prefix(p: str) -> str:
    try:
        p = p.split("/")[-2]
    except IndexError:
        pass

    p = p.replace("work", "")
    p = p.strip("_")

    return p


def get_label(u: mda.Universe) -> str:
    A = clean_universe_prefix(u.filename)

    if False:
        t = u.trajectory.totaltime
        T = "t=%.2fns" % (t / 1000)
        return A + " " + T

    else:
        return process_simulation_name(A)


def process_simulation_name(name: str) -> str:
    """
    Convert internal simulation codes into
    readable names for labels etc...
    """

    number = re.findall(r"\d+-{0,1}\d*", name)

    if number:
        number = number[0]
        if number == "0":
            return "Original"

        return f"Variação #{number}"

    return name


def extract_positions(u: mda.Universe, sel=STANDARD_SELECTION):
    return u.trajectory.timeseries(asel=u.select_atoms(sel))


def pairwise_rmsds(universes: List[mda.Universe]):
    """Compute pairwise RMSDs for position snapshots."""
    size = len(universes)
    rmsd_matrix = np.zeros(shape=(size, size))
    for i, universe_i in enumerate(universes):
        for j, universe_j in enumerate(universes):
            if i == j:
                rmsd = 0
            else:
                align.AlignTraj(
                    universe_i,
                    universe_i,
                    select=STANDARD_SELECTION,
                    in_memory=True
                )

                align.AlignTraj(
                    universe_i,
                    universe_j,
                    select=STANDARD_SELECTION,
                    in_memory=True
                )

                pos_i = extract_positions(universe_i)
                pos_j = extract_positions(universe_j)

                rmsd = rms.rmsd(
                    pos_i.mean(axis=1)[:, None, :],
                    pos_j.mean(axis=1)[:, None, :],
                )

                # rmsd = rms.rmsd(
                #     pos_i,
                #     pos_j
                # )

            rmsd_matrix[i, j] = rmsd
            rmsd_matrix[j, i] = rmsd

    return rmsd_matrix


def pairwise_rmsds_traj(universes: List[mda.Universe], labels: List[str]):
    """Calculate pairwise RMSD for multiple trajectories using Hausdorff distance."""

    # Align trajectories because there seems to be a bug
    # in PSAnalysis.
    for k in universes[1:]:
        align.AlignTraj(
            k,
            universes[0],
            select=STANDARD_SELECTION,
            in_memory=True
        )

    ps = psa.PSAnalysis(universes,
                        labels=labels,
                        reference=universes[0],
                        ref_frame=0,
                        select=STANDARD_SELECTION,
                        path_select=STANDARD_SELECTION)

    ps.generate_paths(align=False, save=False, weights='mass')
    ps.run(metric='hausdorff')
    return ps.D


def load_simulation_prefixes(arguments):

    simulation_prefixes = []
    for autodetect_dir in arguments.AutoDetect:
        simulation_prefixes += autodetect_files(autodetect_dir)

    if arguments.FilePrefix is not None:
        simulation_prefixes += arguments.FilePrefix

    if not simulation_prefixes:
        print("FATAL: No prefixes found.")
        sys.exit(1)

    return simulation_prefixes


def ask_simulation_prefixes(simulation_prefixes):
    print("File prefixes found:")
    for i, prefix in enumerate(simulation_prefixes):
        print(f"{i + 1}:\t" + prefix)

    print("Select all prefixes? (empty input)")
    print("Or input comma separated numbers and " +
          "dash separated intervals to select prefixes.")

    return input(">")


def select_simulation_prefixes(simulation_prefixes, input_string):
    range_descriptors = [
        v.strip()
        for v in input_string.split(",")
    ]

    OutputPrefixes = []

    if not range_descriptors:
        return simulation_prefixes

    for v in range_descriptors:
        try:
            if "-" in v:
                limits = v.split("-")
                assert len(limits) == 2
                F, T = [int(k) for k in limits]
                V = list(range(F, T + 1))
            else:
                V = [int(v)]

        except (ValueError, AssertionError):
            raise Exception("Invalid input.")

        for prefix_idx in V:
            OutputPrefixes.append(simulation_prefixes[prefix_idx - 1])

    for prefix in OutputPrefixes:
        print('\t' + prefix)

    return OutputPrefixes


def concat_filepath(specifiers: List[str]) -> str:
    return "_".join(specifiers) + ".png"


def build_filepath(
        base: str,
        specifiers: List[str],
        arguments) -> Optional[str]:

    if arguments.WriteOutput:
        if arguments.AutoDetect:
            base = arguments.AutoDetect[0]
            if arguments.FilePrefix:
                base += "+"
        else:
            base = "analysis"

        ID = arguments.OutputIdentifier
        if ID:
            specifiers.append(ID)

        return concat_filepath(specifiers)

    return None


def load_universe(simulation_prefix, traj_suffix):

    universe = mda.Universe(
        simulation_prefix + ".gro",
        simulation_prefix + traj_suffix + ".trr"
    )

    align_traj(universe)
    return universe


def align_traj(universe):
    align.AlignTraj(
        universe,
        universe,
        select=STANDARD_SELECTION,
        in_memory=True
    ).run()


def show_universe_information(U: mda.Universe):
    """Show some basic stats for an Universe."""
    print(f"# Atoms:    {len(U.atoms)}")
    print(f"# Frames:   {len(U.trajectory)}")
    print(f"# Residues: {len(U.residues)}")


def global_analysis(arguments):
    """Main pipeline for analysis methods on a trajectory group."""

    simulation_prefixes = load_simulation_prefixes(arguments)

    if arguments.SimulationSelection:
        user_input = arguments.SimulationSelection
    else:
        user_input = ask_simulation_prefixes(simulation_prefixes)

    simulation_prefixes = select_simulation_prefixes(
        simulation_prefixes,
        user_input
    )

    operation_mode = OperationMode(arguments)
    base_filepath = arguments.WriteOutput if arguments.WriteOutput else ""

    print("Data loading done.")

    print("Processing timeseries RMSD plots.")
    session = AnalysisSession()
    for i, simulation_prefix in enumerate(simulation_prefixes):
        print(f"Processsing {i + 1} of {len(simulation_prefixes)}: {simulation_prefix}")
        universe = load_universe(simulation_prefix, arguments.TrajSuffix)

        show_universe_information(universe)

        session.labels.append(get_label(universe))
        if operation_mode.compare_timeseries:
            session.rmsd_series.append(time_series_rmsd(universe, arguments))
            session.rmsf_series.append(time_series_rmsf(universe))

            session.pca_series.append(analyze_pca(universe))

            session.sasa.append(analyze_sasa(universe))

        # Store total time in nanoseconds;
        session.total_times.append(universe.trajectory.totaltime / 1000)

        if operation_mode.compare_pairwise:
            sample = extract_slice_representation(universe)
            session.samples.append(sample)
            session.sample_rmsf.append(time_series_rmsf(sample))

        universe.trajectory.close()
        del universe

        if operation_mode.compare_timeseries:
            plot_series(arguments, base_filepath, session, SeriesMode.STACKED)
            plot_series(arguments, base_filepath, session, SeriesMode.MONOLITHIC)

        if operation_mode.compare_pairwise:
            plot_rmsd_matrices(arguments, base_filepath, session)


def plot_rmsd_matrices(arguments, base_filepath, session):
    selection_frames = list(map(extract_slice_representation, session.samples))
    rmsd_matrix = pairwise_rmsds(selection_frames)

    print(rmsd_matrix)
    mdplots.show_matrix(
        rmsd_matrix,
        session.labels,
        build_filepath(base_filepath, ["pairwise", "rmsds"], arguments)
    )

    rmsd_matrix_traj = pairwise_rmsds_traj(session.samples, session.labels)
    mdplots.show_matrix(
        rmsd_matrix_traj,
        session.labels,
        build_filepath(base_filepath, ["pairwise", "rmsds", "traj"], arguments)
    )


def plot_series(arguments, base_filepath, session, series_mode=SeriesMode.MONOLITHIC):
    series_qualifiers = {
        SeriesMode.MONOLITHIC: {
            "name_appendix": ["mono"],
            "function": mdplots.show_rms_series_monolithic
        },
        SeriesMode.STACKED: {
            "name_appendix": [],
            "function": mdplots.show_rms_series_stacked
        }
    }

    Q = series_qualifiers[series_mode]

    def plot(Q, data, name_segments, mode):
        Q["function"](
            data,
            session.labels,
            session.total_times,
            build_filepath(base_filepath, name_segments + Q["name_appendix"], arguments),
            mode
        )

    plot(Q, session.sample_rmsf, ["ts", "rmsf", "short_sample"], "RMSF")
    plot(Q, session.rmsd_series, ["ts", "rmsd"], "RMSDt")
    plot(Q, session.rmsf_series, ["ts", "rmsf"], "RMSF")
    plot(Q, session.pca_series, ["ts", "variance"], "PCA")
    plot(Q, session.sasa, ["ts", "sasa"], "SASA")


def extract_slice_representation(
        u: mda.Universe,
        slice_position: Optional[Tuple[int, int]] = None
) -> mda.Universe:
    trajectory_length = len(u.trajectory)

    def convert_to_frame(pct, L):
        return int(pct * L)

    slice_position = (
        convert_to_frame(0.85, trajectory_length),
        convert_to_frame(0.9, trajectory_length)
    )

    new_u = mda.Universe(u.filename, u.trajectory.filename)
    new_u.transfer_to_memory(
        start=slice_position[0],
        stop=slice_position[1],
        verbose=False
    )
    align_traj(new_u)
    return new_u


def snapshot_to_universe(source_universe, new_trajectory) -> mda.Universe:
    """Converts a frame snapshot into a universe."""
    # Make a reference structure (need to reshape into a
    # 1-frame "trajectory").

    universe = source_universe.load_new(new_trajectory)

    #universe.trajectory.filename = source_universe.trajectory.filename
    return universe


def align_universe(u: mda.Universe,
                   align_type: AlignType,
                   sel: str = STANDARD_SELECTION):

    align.AlignTraj(
        u,
        u,
        select=sel,
        in_memory=True
    ).run()

    return u


def time_series_rmsd(universe: mda.Universe, arguments, verbose=False) -> List[float]:
    """Extracts the timeseries RMSD from a Universe."""
    rmsds = []

    J = len(universe.trajectory)
    atoms = universe.select_atoms(STANDARD_SELECTION)

    u = align_universe(universe, AlignType.FIRST_FRAME)
    ref = AlignType.FIRST_FRAME.extract(extract_positions(u)).reshape(-1, 3)
    for t, traj in enumerate(u.trajectory):
        if verbose:
            print(f"{t} of {J}")

        frame_rmsd = rms.rmsd(ref, atoms.positions)
        rmsds.append(frame_rmsd)

    # os.remove("rmsfit.xtc")
    return rmsds


def time_series_rmsf(u, sel=STANDARD_SELECTION) -> List[float]:
    """Extract RMSF timeseries from an Universe."""

    align_universe(u, AlignType.MEAN_FRAME)

    rmsf = rms.RMSF(u.select_atoms(sel=sel)).run().rmsf

    return cast(List[float], rmsf)


def analyze_pca(u: mda.Universe, n_dimensions=40):
    """Fetch PCA component contribution values for a single trajectory."""
    PCA = pca.PCA(u, select='backbone')
    space = PCA.run()

    space_3 = space.transform(u.select_atoms('backbone'), 3)
    w = pca.cosine_content(space_3, 0)
    print(w)

    return [
        space.variance[:n_dimensions],
        space.cumulated_variance[:n_dimensions]
    ]


def get_radius(atom):
    radii = {
        "H": 1.1,  # Hydrogen
        "N": 1.6,  # Nitrogen
        "C": 1.7,  # Carbon
        "O": 1.4,  # Oxygen
        "S": 1.8   # Sulfur
    }

    for symbol, radius in radii.items():
        if symbol in atom.name:
            return radius

    return 0


def analyze_sasa(u: mda.Universe):
    atoms = u.select_atoms(STANDARD_SELECTION)
    positions = u.trajectory.timeseries(asel=atoms)

    trajectory_sasa = []
    atom_radius = list(map(get_radius, atoms))
    for frame in np.swapaxes(positions, 0, 1):
        sasa = freesasa.calcCoord(frame.reshape(-1), atom_radius).totalArea()
        trajectory_sasa.append(sasa)

    return np.array(trajectory_sasa)


def main():
    """Executable entrypoint."""
    arguments = parse_arguments()
    global_analysis(arguments)


if __name__ == "__main__":
    main()
