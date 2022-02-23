from typing import List, Optional, cast, Tuple
import enum
import copy
import argparse
import sys
import os
import re
import warnings

import numpy as np

import MDAnalysis as mda
from MDAnalysis.analysis import align, rms, pca

from . import mdplots

warnings.filterwarnings("ignore", category=DeprecationWarning)


def parse_arguments():
    """Parse commandline arguments."""
    parser = argparse.ArgumentParser()

    parser.add_argument("-f", dest='FilePrefix', nargs="*")
    parser.add_argument("-d", dest='AutoDetect', nargs="*")

    parser.add_argument("-t", dest='TrajSuffix', default="")
    parser.add_argument("-M", dest='DoMatrix', action="store_true")
    parser.add_argument("-T", dest='DoTimeseries', action="store_true")

    parser.add_argument("-w", dest='WriteOutput', action="store_true")
    parser.add_argument("-i", "--identifier",
                        dest='OutputIdentifier', required=True)

    parser.add_argument('-m', dest='ReferenceMean', action="store_true")

    parser.add_argument('-s', dest='SimulationSelection')
    return parser.parse_args()


def autodetect_files(root_path, pattern="md.gro") -> List[str]:
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
    NB = re.findall(r"\d+-{0,1}\d*", name)

    if NB:
        number = NB[0]
        if number == "0":
            return "Original"

        return f"Variação #{number}"

    return name

def RMSDStudy(us, unames):
    POS = []

    traj_idx = [0, -1]
    ATOM_ID = "name CA"
    ATOM_ID = "backbone"

    for i, universe in enumerate(us):
        bb = universe.select_atoms(ATOM_ID)
        for j in traj_idx:
            universe.trajectory[j]
            w = Positions(bb.positions.copy(), unames[i], j)
            POS.append(w)

    return POS


def pairwise_rmsds(POS):
    SIZE = len(POS)
    w = np.zeros(shape=(SIZE, SIZE))
    for i, pi in enumerate(POS):
        for j, pj in enumerate(POS):
            if i == j:
                v = 0
            else:
                v = rms.rmsd(pi.positions, pj.positions)
            w[i, j] = v
            w[j, i] = v

    return w


def loadSimulationPrefixes(arguments):

    simulation_prefixes = []
    for AutoDetectDir in arguments.AutoDetect:
        simulation_prefixes += autodetect_files(AutoDetectDir)

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

    U = mda.Universe(
        simulation_prefix + ".gro",
        simulation_prefix + traj_suffix + ".trr"
    )

    align_traj(U)
    return U


def align_traj(universe):
    align.AlignTraj(
        universe,
        universe,
        select='name CA',
        in_memory=True
    ).run()


def show_universe_information(U: mda.Universe):
    print(f"# Atoms:  {len(U.atoms)}")
    print(f"# Frames: {len(U.trajectory)}")


def analyzeMD(arguments):

    simulation_prefixes = loadSimulationPrefixes(arguments)

    if arguments.SimulationSelection:
        user_input = arguments.SimulationSelection
    else:
        user_input = ask_simulation_prefixes(simulation_prefixes)

    simulation_prefixes = select_simulation_prefixes(
        simulation_prefixes,
        user_input
    )

    # us = map(lambda sp: load_universe(sp, arguments), simulation_prefixes)
    # can access via segid (4AKE) and atom name
    # we take the first atom named N and the last atom named C
    # nterm = u.select_atoms('segid 4AKE and name N')[0]
    # cterm = u.select_atoms('segid 4AKE and name C')[-1]

    base_filepath = arguments.WriteOutput if arguments.WriteOutput else None

    print("Data loading done.")

    if arguments.DoTimeseries:
        print("Processing timeseries RMSD plots.")
        labels = []
        rmsd_series = []
        rmsf_series = []
        pca_series = []
        total_times = []

        samples = []

        for i, SP in enumerate(simulation_prefixes):
            print(f"Processsing {i + 1} of {len(simulation_prefixes)}: {SP}")
            u = load_universe(SP, arguments.TrajSuffix)

            show_universe_information(u)

            labels.append(get_label(u))
            rmsd_series.append(time_series_rmsd(u, arguments))
            rmsf_series.append(time_series_rmsf(u))

            pca_series.append(analyze_pca(u))

            # Store total time in nanoseconds;
            total_times.append(u.trajectory.totaltime / 1000)

            #samples.append(extract_slice_representation(u))
            u.trajectory.close()
            del u

        if True:
            mdplots.show_rms_series_stacked(
                rmsd_series,
                labels,
                total_times,
                build_filepath(base_filepath, ["ts", "rmsd"], arguments),
                "RMSDt"
            )

            mdplots.show_rms_series_stacked(
                rmsf_series,
                labels,
                total_times,
                build_filepath(base_filepath, ["ts", "rmsf"], arguments),
                "RMSF"
            )

            mdplots.show_rms_series_stacked(
                pca_series,
                labels,
                total_times,
                build_filepath(base_filepath, ["ts", "variance"], arguments),
                "PCA"
            )


        if False:
            mdplots.show_rms_series_monolithic(
                rmsd_series,
                labels,
                total_times,
                build_filepath(base_filepath, ["tsmono", "rmsd"], arguments),
                "RMSDt"
            )

            mdplots.show_rms_series_monolithic(
                rmsf_series,
                labels,
                total_times,
                build_filepath(base_filepath, ["tsmono", "rmsf"], arguments),
                "RMSF"
            )


def extract_slice_representation(u: mda.Universe, slice_position: Optional[Tuple[int, int]] = None):
    L = len(u._trajectory)

    def convert_to_frame(pct, L):
        return int(pct * L)

    slice_position = (
        convert_to_frame(0.85, L),
        convert_to_frame(0.9, L)
    )

    new_u = snapshot_to_universe(
        u.atoms,
        [k for k in u.trajectory[slice(*slice_position)]]
    )

    align_traj(new_u)
    AlignType.MEAN_FRAME.extract(
        new_u.trajectory.timeseries(asel="protein and name CA")
    )
    return new_u


class AlignType(enum.Enum):
    """
    Extracts a single frame structural summary from
    trajectories by using different methods.
    """
    FIRST_FRAME = 0
    MEAN_FRAME = 1

    def extract(self, traj):
        if self == AlignType.FIRST_FRAME:
            return traj[:, 0, :][:, None, :]
        if self == AlignType.MEAN_FRAME:
            return traj.mean(axis=1)[:, None, :]


def snapshot_to_universe(atoms, snapshot):
    """Converts"""
    # Make a reference structure (need to reshape into a
    # 1-frame "trajectory").
    return mda.Merge(atoms).load_new(
        snapshot,
        order="afc"
    )


def align_universe(u: mda.Universe,
                   align_type: AlignType,
                   sel: str = "protein and name CA"):

    atoms = u.select_atoms(sel)

    ref_coordinates = align_type.extract(u.trajectory.timeseries(asel=atoms))

    ref = snapshot_to_universe(atoms, ref_coordinates)

    align.AlignTraj(
        u,
        ref,
        select=sel,
        in_memory=True
    ).run()

    return ref, u


def time_series_rmsd(u, arguments, verbose=False) -> List[float]:
    rmsds = []

    J = len(u.trajectory)
    atoms = u.select_atoms("protein and name CA")

    ref, u = align_universe(u, AlignType.FIRST_FRAME)
    ref_atoms = ref.select_atoms("protein and name CA")
    for t, traj in enumerate(u.trajectory):
        if verbose:
            print(f"{t} of {J}")

        frame_rmsd = rms.rmsd(ref_atoms.positions, atoms.positions)
        rmsds.append(frame_rmsd)

    # os.remove("rmsfit.xtc")

    return rmsds


def time_series_rmsf(u, end_pct=100, sel="protein and name CA") -> List[float]:

    _, aligned_u = align_universe(u, AlignType.MEAN_FRAME)

    rmsf = rms.RMSF(aligned_u.select_atoms(sel=sel)).run().rmsf

    return cast(List[float], rmsf)


def analyze_pca(u: mda.Universe, n=40):
    PCA = pca.PCA(u, select='backbone')
    space = PCA.run()

    space_3 = space.transform(u.select_atoms('backbone'), 3)
    w = pca.cosine_content(space_3, 0)
    print(w)

    return [
        space.variance[:n],
        space.cumulated_variance[:n]
    ]


def main():
    arguments = parse_arguments()
    analyzeMD(arguments)


if __name__ == "__main__":
    main()
