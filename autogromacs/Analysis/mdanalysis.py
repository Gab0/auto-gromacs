from typing import List, Optional, cast, Tuple

import enum
import sys
import os
import warnings
import pickle

import Bio.PDB.Polypeptide as polyp
import freesasa
import numpy as np

from sklearn.preprocessing import normalize

import MDAnalysis as mda
from MDAnalysis.analysis import align, rms, pca, psa

from MDAnalysis.analysis.dihedrals import Ramachandran

from . import cli_arguments, mdplots, user_input
from . import crosscorr, dimension_reduction, superposition
from . import mutation_labels

from antigen_protocol.Mutation import structure_name
from antigen_protocol.ProteinSequence import Antigens
from antigen_protocol.StructureUtils import BasicStructureOperations as BSO

warnings.filterwarnings("ignore", category=DeprecationWarning)

STANDARD_SELECTION = "protein and name CA"


class Positions():
    def __init__(self, pos, universe_name, traj_idx):
        self.positions = pos
        self.universe_name = universe_name
        self.traj_idx = traj_idx

        self.label = self.make_label()

    def make_label(self):
        return "_".join([
            mutation_labels.clean_universe_prefix(self.universe_name),
            index_label(self.traj_idx)
        ])


class AnalysisSession():
    """Holds what is saved from all loaded Universes for a single analysis."""
    store_universe: bool = False
    plot_suffix: List[str]
    sample_pct: Optional[Tuple[float, float]] = None

    universes: List[mda.Universe]
    labels: List[str]
    rmsd_series: List[List[float]]
    rmsf_series: List[List[float]]
    pca_series: List[np.ndarray]
    total_times: List[int]
    sasa: List[np.ndarray]
    radgyr: List[List[float]]
    snapshots: List[np.ndarray]
    secondary_structure_n: List[List[float]]
    sequences: List[str]

    def __init__(self, store_universe, plot_suffix, sample_pct, selector=None):
        self.store_universe = store_universe
        self.plot_suffix = plot_suffix
        self.sample_pct = sample_pct
        self.selector = selector
        self.universes = []

        self.secondary_structure_n = []

        for feature_name, _ in self.features:
            self.__dict__[feature_name] = []

    @property
    def features(self):
        """Declare all the features that can be extracted from"""
        return [
            ("labels", mutation_labels.get_label),
            # Total times in nanoseconds;
            ("total_times", lambda u: u.trajectory.totaltime / 1000),
            ("rmsd_series", time_series_rmsd),
            ("rmsf_series", time_series_rmsf),
            #("pca_series", analyze_pca),
            ("sasa", analyze_sasa),
            ("radgyr", analyze_radgyr),
            ("secondary_structure_n", None),
            ("sequences", structure_sequence)
        ]

    def update(self, universe: mda.Universe, simulation_directory):
        """Add analysis for a single Universe into this session."""

        print(f"Updating. {self.store_universe}")

        if self.sample_pct is not None:
            universe = extract_slice_representation(universe, self.sample_pct)
        elif self.sample_pct == (1.0, 1.0):
            universe = get_best_stable_window(universe)

        if self.store_universe:
            self.universes.append(universe)

        selector = STANDARD_SELECTION
        if self.selector is not None:
            selector += " and " + self.selector

        for feature_name, feature_extractor in self.features:
            if feature_extractor is not None:
                if "selector" in feature_extractor.__code__.co_varnames:
                    feature = feature_extractor(universe, selector)
                else:
                    feature = feature_extractor(universe)
                self.__dict__[feature_name].append(feature)

        secondary_n = analyze_secondary(simulation_directory)

        print("SEC")
        print(secondary_n)

        (sf, st) = slice_to_indexes(self.sample_pct, len(secondary_n))

        secondary_n_output = secondary_n[sf:st]
        print("N OUT")
        print(secondary_n_output)
        self.secondary_structure_n.append(secondary_n_output)
        #self.snapshots.append()

    def check(self):
        """Some checks for session integrity."""
        for key, value in self.__dict__.items():
            if isinstance(value, list):
                print(len(value))

    def check_trajectory_stability(self) -> bool:
        """Checks if all stored trajectory fragments are stable."""
        stable = True
        for label, rmsd_traj in zip(self.labels, self.rmsd_series):
            k = [max(rmsd_traj), min(rmsd_traj)]
            if abs(k[0] - k[1]) > 5:
                print(f"Unstable trajectory for {label}.")
                stable = False

        return stable

    def select_simulation_indexes(self, selected_indexes: List[int]) -> None:
        """
        Select the subset of stored data that
        corresponds to a given simulation index set.
        """
        for feature_name, _ in self.features:
            self.__dict__[feature_name] = [
                self.__dict__[feature_name][idx]
                for idx in selected_indexes
            ]


class SeriesMode(enum.Enum):
    """Different types of series visualisation."""
    MONOLITHIC = 0
    STACKED = 1


class OperationMode():
    """
    Handles the operation mode for the anaylsis,
    which is based on the user input arguments.
    """
    compare_pairwise = True
    compare_full_timeseries = True
    compare_short_timeseries = False
    compare_samples = True

    def __init__(self, arguments):
        self.compare_full_timeseries = not arguments.matrix_only


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

        return traj


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


def index_label(label_id: int) -> str:
    """"""
    common = {
        0: "initial",
        -1: "last"
    }
    try:
        return common[label_id]
    except KeyError:
        return str(label_id)


def extract_positions(universe: mda.Universe, sel=STANDARD_SELECTION):
    """ Extract trajectory atomic position vectors from a Universe."""
    return universe.trajectory.timeseries(asel=universe.select_atoms(sel))


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


def normalize_rmsf(rmsf_series):
    return normalize(rmsf_series, axis=1, norm='max')


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


def concat_filepath(specifiers: List[str]) -> str:
    return "_".join(specifiers) + ".png"


def build_filepath(
        base: str,
        specifiers: List[str],
        arguments) -> Optional[str]:

    if arguments.no_plot:
        return None

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


def session_selector(arguments, sessions: List[AnalysisSession]) -> List[AnalysisSession]:
    """Select simulations for an existing session."""

    session = sessions[0]

    selected = arguments.SimulationSelection
    if not selected:
        selected = user_input.ask_simulation_prefixes(session.labels)

    selected_indexes = user_input.process_range_descriptors(
        selected,
        len(session.labels)
    )

    for session in sessions:
        session.select_simulation_indexes(selected_indexes)

    return sessions


def determine_sessions(operation_mode: OperationMode) -> List[AnalysisSession]:
    sessions = [AnalysisSession(True, ["total"], None)]

    if operation_mode.compare_short_timeseries:
        sessions += [
            #AnalysisSession(True, ["short-sample-080"], (0.8, 0.85)),
            #AnalysisSession(True, ["short-sample-085"], (0.85, 0.9)),
            #AnalysisSession(True, ["short-sample-090"], (0.9, 0.95)),
            #AnalysisSession(True, ["short-sample-best"], (1.0, 1.0))
        ]

    if operation_mode.compare_full_timeseries:
        sessions += [
            # RESID 129 divides the two domains in SRS29B;
            AnalysisSession(True, ["total", "A"], None, selector="resid 1:129"),
            AnalysisSession(True, ["total", "B"], None, selector="resid 129:3000")
       ]

    return sessions


def global_analysis(arguments):
    """Main pipeline for analysis methods on a trajectory group."""

    simulation_prefixes = load_simulation_prefixes(arguments)

    if arguments.SimulationSelection:
        user_selection = arguments.SimulationSelection
    else:
        user_selection = user_input.ask_simulation_prefixes(simulation_prefixes)

    if user_selection != "":
        simulation_prefixes = user_input.select_simulation_prefixes(
            simulation_prefixes,
            user_selection
        )

    operation_mode = OperationMode(arguments)

    sessions = determine_sessions(operation_mode)

    print("Data loading done.")

    print("Processing timeseries RMSD plots.")

    for i, simulation_prefix in enumerate(simulation_prefixes):
        print(
            f"Processsing {i + 1} of {len(simulation_prefixes)}: "
            + f"{simulation_prefix}"
        )

        universe = load_universe(simulation_prefix, arguments.TrajSuffix)
        simulation_directory = os.path.join(*os.path.split(simulation_prefix)[:-1])
        show_universe_information(universe)

        for session in sessions:
            session.update(universe, simulation_directory)
            if session.universes:
                print("Building snapshot...")
                superposition.build_snapshot(
                    session.universes,
                    [arguments.OutputIdentifier] + session.plot_suffix,
                    session.labels
                )
                for _universe in session.universes:
                    _universe.trajectory.close()
                    del _universe
            else:
                print("Session has no attached universes.")

        # Close full-length universe to preserve RAM memory.
        universe.trajectory.close()
        del universe

    for session in sessions:
        # FIXME: Is `del` required here?
        del session.universes
        session.universes = []

    return sessions


def plot_sessions(sessions, arguments):
    """Plot routines for all gathered data."""

    operation_mode = OperationMode(arguments)

    # TODO: Source this from an argument?
    base_filepath = ""

    for session in sessions:
        session.check_trajectory_stability()

        plot_series(
            arguments,
            base_filepath,
            session,
            session.plot_suffix,
            series_mode=SeriesMode.STACKED
        )

        plot_series(
            arguments,
            base_filepath,
            session,
            session.plot_suffix,
            series_mode=SeriesMode.MONOLITHIC
        )

        if operation_mode.compare_pairwise:
            if session.universes:
                plot_rmsd_matrices(arguments, base_filepath, session)

        for universe, label in zip(session.universes, session.labels):
            crosscorr.trajectory_cross_correlation(
                universe,
                build_filepath(base_filepath, ["cross-corr", label], arguments)
            )

            ramachandran(
                universe,
                label,
                build_filepath(base_filepath, ["rama", label], arguments)
            )


def plot_umap(session: AnalysisSession, arguments, base_filepath: str):
    try:
        rmsf_norm = normalize_rmsf(session.rmsf_series)

        rmsf_2d = dimension_reduction.umap_reduce(rmsf_norm)
        dimension_reduction.plot_2D(
            rmsf_2d,
            session.labels,
            build_filepath(base_filepath, ["umap-rmsf"] + session.plot_suffix, arguments)
        )

        rmsd_2d = dimension_reduction.umap_reduce(session.rmsd_series)
        dimension_reduction.plot_2D(
            rmsd_2d,
            session.labels,
            build_filepath(base_filepath, ["umap-rmsd"] + session.plot_suffix, arguments)
        )
    except ValueError as e:
        print("UMAP failure!")
        print(e)


def ramachandran(u: mda.Universe, label: str, output_filepath: str):

    sel = u.select_atoms(STANDARD_SELECTION)
    R = Ramachandran(sel).run()
    mdplots.plot_ramachandran(
        R,
        label,
        output_filepath
    )


def plot_rmsd_matrices(arguments, base_filepath, session):

    rmsd_matrix = pairwise_rmsds(session.universes)

    print(rmsd_matrix)
    mdplots.show_matrix(
        rmsd_matrix,
        session.labels,
        build_filepath(
            base_filepath,
            ["pairwise", "rmsds"] + session.plot_suffix,
            arguments
        )
    )

    rmsd_matrix_traj = pairwise_rmsds_traj(session.universes, session.labels)
    mdplots.show_matrix(
        rmsd_matrix_traj,
        session.labels,
        build_filepath(
            base_filepath,
            ["pairwise", "rmsds", "traj"] + session.plot_suffix,
            arguments
        )
    )


def load_mutation_labels(session, reference_structure_path: Optional[str]) -> Optional[List[str]]:

    if reference_structure_path is None:
        return None

    reference_sequence = BSO.read_structure_sequence(reference_structure_path)

    print("Loaded reference sequence:")
    print(reference_sequence)
    # assert session.sequences[0] != session.sequences[1]

    mutation_vectors = [
        Antigens.CreateMutationVector(
            [seq],
            reference_sequence,
            [label]
        )
        for seq, label in zip(session.sequences, session.labels)
    ]

    for mv in mutation_vectors:
        print(f">{[m for m in mv if m]}")

    session_mutations = [
        [mut for mut in mutation_vector if mut]
        for mutation_vector in mutation_vectors
    ]

    extra_labels = [
        " ".join([mut.show_mutations()[0] for mut in session_muts])
        for session_muts in session_mutations
    ]

    return extra_labels


def plot_series(
        arguments,
        base_filepath,
        session,
        extra_identifier=[],
        series_mode=SeriesMode.MONOLITHIC):

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

    # Create extra labels for all structures.
    extra_labels = load_mutation_labels(session, arguments.reference_structure)


    def plot(Q, data, name_segments, mode):
        try:
            Q["function"](
                data,
                session.labels,
                session.total_times,
                build_filepath(
                    base_filepath,
                    name_segments + Q["name_appendix"] + extra_identifier,
                    arguments),
                mode,
                extra_labels
            )
        except Exception as exception:
            print(f"Could not create {mode} plots for {extra_identifier}")
            print(f"Reason: {exception}\n")

    # for feature, _, plot in session.plotable_features:
    #     plot(Q, session.getattr(feature),)
    plot(Q, session.rmsd_series, ["ts", "rmsd"], "RMSDt")
    plot(Q, session.rmsf_series, ["ts", "rmsf"], "RMSF")
    #plot(Q, session.pca_series, ["ts", "variance"], "PCA")
    plot(Q, session.sasa, ["ts", "sasa"], "SASA")
    plot(Q, session.radgyr, ["ts", "radgyr"], "RADGYR")

    # Process and plot 'secondary struct n' data structures.
    sec_struct = session.secondary_structure_n
    # print(f"Secondary structure data shape: {identify_list_object(sec_struct)}")
    plot(Q, analyze_secondary_structs(sec_struct), ["ts", "secondary", "strut"], "NSECONDARY")


def analyze_secondary_structs(obj):
    """
    Process 'secondary struct n' data.
    Raw, primarily parsed data from source
    (GROMACS's do_dssp's xcount.xvg file) is stored in the
    sessions in order to maximize compatibility with different GROMACS versions.
    A downside further processing is required to plot the information,
    and this is done here.

    The post-processing consists in simply summing the different
    columns in the original data.
    But this may require updates in the future.
    """

    return [
        [sum(j[1:]) for j in k.T]
        for k in obj
    ]


def identify_list_object(obj):
    """ Identify shapes on possibly 'ragged list' arrangements."""

    try:
        return obj.shape
    except AttributeError:
        return [k.shape for k in obj]


def slice_to_indexes(
        slice_pct: Optional[Tuple[float, float]],
        total_length: int
) -> Tuple[int, int]:
    if slice_pct is None:
        return (0, total_length)

    def convert_to_index(pct, L):
        return int(pct * L)

    pct_start, pct_end = slice_pct
    return (
        convert_to_index(pct_start, total_length),
        convert_to_index(pct_end, total_length)
    )


def extract_slice_representation(
        u: mda.Universe,
        slice_position_pct: Tuple[float, float] = (0.85, 0.9)
) -> mda.Universe:
    trajectory_length = len(u.trajectory)

    slice_position = slice_to_indexes(slice_position_pct, trajectory_length)

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

    # universe.trajectory.filename = source_universe.trajectory.filename
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


def structure_sequence(universe: mda.Universe) -> str:
    protein_string = ""
    for res in universe.residues:
        try:
            protein_string += polyp.three_to_one(res.resname)
        except KeyError:
            continue

    return protein_string


def time_series_rmsd(universe: mda.Universe, selector: str, verbose=False) -> List[float]:
    """Extracts the RMSD timeseries from a Universe."""

    rmsds = []

    total_length = len(universe.trajectory)

    aligned_universe = align_universe(universe, AlignType.FIRST_FRAME, sel=selector)
    atoms = aligned_universe.select_atoms(selector)
    ref = AlignType.FIRST_FRAME.extract(extract_positions(aligned_universe, sel=selector)).reshape(-1, 3)
    for frame_idx, _ in enumerate(aligned_universe.trajectory):
        if verbose:
            print(f"{frame_idx} of {total_length}")
            if not frame_idx:
                print("RMSD shapes:")
                print(ref.shape)
                print(atoms.positions.shape)

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
    pca_analysis = pca.PCA(u, select='backbone')
    space = pca_analysis.run()

    space_3 = space.transform(u.select_atoms('backbone'), 3)
    w = pca.cosine_content(space_3, 0)
    print(w)

    return [
        space.variance[:n_dimensions],
        space.cumulated_variance[:n_dimensions]
    ]


def get_atom_radius(atom):
    """Get atom radii."""
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

    raise Exception("Unknown ATOM {atom.name}.")


def analyze_sasa(u: mda.Universe) -> np.ndarray:
    """Extract SASA value for each trajectory frame."""
    atoms = u.select_atoms(STANDARD_SELECTION)
    positions = u.trajectory.timeseries(asel=atoms)

    trajectory_sasa = []
    atom_radius = list(map(get_atom_radius, atoms))
    for frame in np.swapaxes(positions, 0, 1):
        sasa = freesasa.calcCoord(frame.reshape(-1), atom_radius).totalArea()
        trajectory_sasa.append(sasa)

    return np.array(trajectory_sasa)


def analyze_radgyr(u: mda.Universe) -> List[float]:
    """Extract the radius of gyration metric for each trajectory frame."""
    trajectory_radgyr = []
    atoms = u.select_atoms(STANDARD_SELECTION)
    for _ in u.trajectory:
        trajectory_radgyr.append(atoms.radius_of_gyration())

    return trajectory_radgyr


def analyze_secondary(simulation_directory: str) -> List[float]:
    """
    Extract the number of residues participating in secondary structures
    on each frame.
    """
    xvg = os.path.join(simulation_directory, "scount.xvg")
    structs = np.loadtxt(xvg, comments=["@", "#"], unpack=True)

    if isinstance(structs[1], List):
        structural = structs[1]
    else:
        structural = structs

    return cast(List[float], structural)


def get_best_stable_window(
        universe: mda.Universe,
        starting_point: float = 0.7,
        window_size: float = 0.05) -> mda.Universe:
    """Get the most stable window based on RMSD from an Universe."""
    window_frames = ((s, s + window_size) for s in np.arange(starting_point, 1 - window_size, 0.1))
    frames = []
    delta_rmsds = []

    def evaluate_rmsd(rmsd):
        return np.abs(np.min(rmsd) - np.max(rmsd))

    for frame_bounds in window_frames:
        frame = extract_slice_representation(universe, frame_bounds)
        frames.append(frame)
        delta_rmsds.append(evaluate_rmsd(time_series_rmsd(universe, frame)))

    return frames[delta_rmsds.index(min(delta_rmsds))]


def load_session(session_file) -> List[AnalysisSession]:
    """Load a previously saved session, which is a pickled file."""
    with open(session_file, 'rb') as fin:
        return cast(List[AnalysisSession], pickle.load(fin))


def main():
    """Executable entrypoint."""
    arguments = cli_arguments.parse_arguments()

    if arguments.load_session:
        sessions = load_session(arguments.load_session)
    else:
        sessions = global_analysis(arguments)

    if arguments.write_session:
        with open(arguments.write_session, 'wb') as fout:
            pickle.dump(sessions, fout)

    sessions = session_selector(arguments, sessions)
    if not arguments.no_plot:
        plot_sessions(sessions, arguments)


if __name__ == "__main__":
    main()
