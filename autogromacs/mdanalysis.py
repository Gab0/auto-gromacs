
from typing import List, Optional, cast
import enum
import argparse
import sys
import os
import re
import numpy as np
import warnings

import MDAnalysis as mda
from MDAnalysis.analysis import align, rms, pca

import numpy.linalg

from . import mdplots

warnings.filterwarnings("ignore", category=DeprecationWarning)


def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument("-f", dest='FilePrefix', nargs="*")
    parser.add_argument("-d", dest='AutoDetect', nargs="*")

    parser.add_argument("-t", dest='TrajSuffix', default="")
    parser.add_argument("-M", dest='DoMatrix', action="store_true")
    parser.add_argument("-T", dest='DoTimeseries', action="store_true")
    parser.add_argument("-w", dest='WriteOutput', action="store_true")

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
    NB = re.findall(r"\d+-{0,1}\d*", name)[0]
    if NB == "0":
        return "Original"

    return f"Variação #{NB}"


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

    SimulationPrefixes = []
    for AutoDetectDir in arguments.AutoDetect:
        SimulationPrefixes += autodetect_files(AutoDetectDir)

    if arguments.FilePrefix is not None:
        SimulationPrefixes += arguments.FilePrefix

    if not SimulationPrefixes:
        print("FATAL: No prefixes found.")
        sys.exit(1)

    return SimulationPrefixes


def ask_simulation_prefixes(SimulationPrefixes):
    print("File prefixes found:")
    for i, prefix in enumerate(SimulationPrefixes):
        print(f"{i + 1}:\t" + prefix)

    print("Select all prefixes? (empty input)")
    print("Or input comma separated numbers and " +
          "dash separated intervals to select prefixes.")

    return input(">")


def select_simulation_prefixes(SimulationPrefixes, input_string):
    q = [
        v.strip()
        for v in input_string.split(",")
    ]

    OutputPrefixes = []

    if not q:
        return SimulationPrefixes

    for v in q:
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
            OutputPrefixes.append(SimulationPrefixes[prefix_idx - 1])

    for prefix in OutputPrefixes:
        print('\t' + prefix)

    return OutputPrefixes


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
        return base + "_" + "_".join(specifiers) + ".png"

    return None


def load_universe(SimulationPrefix, arguments):

    U = mda.Universe(
        SimulationPrefix + ".gro",
        SimulationPrefix + arguments.TrajSuffix + ".trr"
    )

    align.AlignTraj(
        U,
        U,
        select='name CA',
        in_memory=True
    ).run()

    return U


def analyzeMD(arguments):

    SimulationPrefixes = loadSimulationPrefixes(arguments)

    if arguments.SimulationSelection:
        user_input = arguments.SimulationSelection
    else:
        user_input = ask_simulation_prefixes(SimulationPrefixes)

    SimulationPrefixes = select_simulation_prefixes(
        SimulationPrefixes,
        user_input
    )

    # us = map(lambda sp: load_universe(sp, arguments), SimulationPrefixes)
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

        for i, SP in enumerate(SimulationPrefixes):
            print(f"Processsing {i + 1} of {len(SimulationPrefixes)}: {SP}")
            u = load_universe(SP, arguments)
            labels.append(get_label(u))
            rmsd_series.append(time_series_rmsd(u, arguments))
            rmsf_series.append(time_series_rmsf(u))

            pca_series.append(analyze_pca(u))

            # Store total time in nanoseconds;
            total_times.append(u.trajectory.totaltime / 1000)

            u.trajectory.close()
            del u

        mdplots.show_rms_series(
            rmsd_series,
            labels,
            total_times,
            build_filepath(base_filepath, ["tsp", "rmsd"], arguments),
            "RMSDt"
        )

        mdplots.show_rms_series_monolithic(
            rmsd_series,
            labels,
            total_times,
            build_filepath(base_filepath, ["tsmono", "rmsd"], arguments),
            "RMSDt"
        )

        mdplots.show_rms_series(
            rmsf_series,
            labels,
            total_times,
            build_filepath(base_filepath, ["ts", "rmsf"], arguments),
            "RMSF"
        )

        mdplots.show_rms_series_monolithic(
            rmsf_series,
            labels,
            total_times,
            build_filepath(base_filepath, ["tsmono", "rmsf"], arguments),
            "RMSF"
        )

        mdplots.show_rms_series_monolithic(
            pca_series,
            labels,
            total_times,
            build_filepath(base_filepath, ["tsmono", "pca"], arguments),
            "PCA"
        )


class AlignType(enum.Enum):
    FIRST_FRAME = 0
    MEAN_FRAME = 1


def align_universe(u: mda.Universe,
                   align_type: AlignType,
                   sel: str = "protein and name CA"):

    atoms = u.select_atoms(sel)

    if align_type == AlignType.MEAN_FRAME:
        ref_coordinates = u.trajectory.timeseries(asel=atoms).mean(axis=1)
    else:
        ref_coordinates = u.trajectory.timeseries(asel=atoms)[:, 0, :]

    REF = ref_coordinates.copy()

    # Make a reference structure (need to reshape into a
    # 1-frame "trajectory").
    ref = mda.Merge(atoms).load_new(
        ref_coordinates[:, None, :],
        order="afc"
    )

    align.AlignTraj(
        u,
        ref,
        select="protein and name CA",
        in_memory=True
    ).run()

    return REF, atoms


def time_series_rmsd(u, arguments, verbose=False) -> List[float]:
    rmsds = []

    J = len(u.trajectory)

    ref, atoms = align_universe(u, AlignType.FIRST_FRAME)
    for t, traj in enumerate(u.trajectory):

        if verbose:
            print(f"{t} of {J}")

        v = rms.rmsd(ref, atoms.positions)
        rmsds.append(v)

    # os.remove("rmsfit.xtc")

    return rmsds


def time_series_rmsf(u, end_pct=100) -> List[float]:

    ref, atoms = align_universe(u, AlignType.MEAN_FRAME)

    rmsf = rms.RMSF(atoms).run().rmsf

    return cast(List[float], rmsf)


def analyze_pca(u: mda.Universe):
    PCA = pca.PCA(u, select='backbone')
    space = PCA.run()

    space_3 = space.transform(u.select_atoms('backbone'), 3)
    w = pca.cosine_content(space_3, 0)
    print(w)

    return space.cumulated_variance[:10]


def main():
    arguments = parse_arguments()
    analyzeMD(arguments)


if __name__ == "__main__":
    main()
