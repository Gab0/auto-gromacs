
from typing import List, Union, Optional
import argparse
import sys
import os
import re
import numpy as np
import warnings
import MDAnalysis
import MDAnalysis as mda
from MDAnalysis.analysis import diffusionmap, align, rms

import numpy.linalg
import seaborn as sns
import matplotlib.pyplot as plt

sns.set_theme(style="darkgrid")
plt.rcParams["axes.labelsize"] = 15
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


def selectSimulationPrefixes(SimulationPrefixes):
    print("File prefixes found:")
    for i, prefix in enumerate(SimulationPrefixes):
        print(f"{i + 1}:\t" + prefix)

    print("Select all prefixes? (empty input)")
    print("Or input comma separated numbers and " +
          "dash separated intervals to select prefixes.")

    q = input(">")

    q = [
        v.strip()
        for v in q.split(",")
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
            print("Invalid input.")
            exit(1)

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
            base = arguments.AutoDetect
            if arguments.FilePrefix:
                base += "+"
        else:
            base = "analysis"
        return base + "_" + "_".join(specifiers) + ".png"

    return None


def load_universe(SimulationPrefix, arguments):

    U = MDAnalysis.Universe(
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
    SimulationPrefixes = selectSimulationPrefixes(SimulationPrefixes)

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

        for i, SP in enumerate(SimulationPrefixes):
            print(f"Processsing {i + 1} of {len(SimulationPrefixes)}: {SP}")
            u = load_universe(SP, arguments)
            labels.append(get_label(u))
            rmsd_series.append(time_series_rmsd(u, arguments))
            rmsf_series.append(time_series_rmsf(u))

            u.trajectory.close()
            del u

        show_rms_series(
            rmsd_series,
            labels,
            build_filepath(base_filepath, ["tsp", "rmsd"], arguments),
            "RMSDt"
        )

        show_rms_series_monolithic(
            rmsd_series, labels,
            build_filepath(base_filepath, ["tsmono", "rmsd"], arguments),
            "RMSDt"
        )

        show_rms_series(
            rmsf_series,
            labels,
            build_filepath(base_filepath, ["ts", "rmsf"], arguments),
            "RMSF"
        )

        show_rms_series_monolithic(
            rmsf_series,
            labels,
            build_filepath(base_filepath, ["tsmono", "rmsf"], arguments),
            "RMSF"
        )


def show_matrix(results, labels, filepath: Union[str, None]):
    fig, ax = plt.subplots()

    im = ax.imshow(results, cmap='viridis')

    fig.colorbar(im, ax=ax, label=r'RMSD ($\AA$)')
    # We want to show all ticks...
    U = np.arange(len(labels))
    ax.set_yticks(U)
    plt.xticks(range(len(results)), labels, rotation='vertical')

    # ax.set_xticks(U, rotation='vertical')
    # ... and label them with the respective list entries
    ax.set_yticklabels(labels)
    ax.set_xticklabels(labels)

    plt.tight_layout()

    if filepath is not None:
        plt.savefig(filepath)
    else:
        plt.show()


def show_rms_series_monolithic(
        rms_series: List[List[float]],
        labels: List[str],
        filepath: Union[str, None],
        mode: str):

    fig, ax = plt.subplots()

    fig.set_figwidth(9.6)

    for i, Xa in enumerate(rms_series):
        ax.plot(range(len(Xa)), Xa)

    YL = r"Distância ($\AA$)"
    if mode == "RMSD":
        XL = "Frame"

    elif mode == "RMSF":
        XL = "Residue"
    else:
        exit(1)

    # ax.set_title(mode)
    ax.set_xlabel(XL)
    ax.set_ylabel(YL)

    ax.legend(labels)
    plt.tight_layout()

    if filepath is not None:
        plt.savefig(filepath)
    else:
        plt.show()


def frames_to_time(frames: Union[List[float], List[int]],
                   total_time: int) -> List[float]:
    total_frames = frames[-1]

    def to_t(v: Union[int, float]) -> float:
        return v / total_frames * total_time

    return list(map(to_t, frames))


def show_rms_series(
        rms_series: List[List[float]],
        labels: List[str],
        filepath: Union[str, None],
        mode: str):

    N = len(labels)
    ncols = 3
    nrows = round(np.ceil(N / ncols))

    assert ncols * nrows >= N
    fig, ax = plt.subplots(nrows, ncols)

    axk = ax.ravel()

    for i, (vals, label) in enumerate(zip(rms_series, labels)):

        Y = vals
        X: Union[List[float], List[int]] = list(range(len(Y)))

        YL = r"Distância ($\AA$)"
        if mode == "RMSDf":
            XL = "Frame"
        elif mode == "RMSDt":
            XL = "Time (ns)"
            X = frames_to_time(X, 64)
        elif mode == "RMSF":
            XL = "Residue"
        else:
            exit(1)

        axk[i].plot(X, Y, "b-")

        axk[i].set_title(label)

        axk[i].set_xlabel(XL)
        axk[i].set_ylabel(YL)

    # plt.title(mode)
    plt.tight_layout()

    if filepath is not None:
        plt.savefig(filepath)
    else:
        plt.show()


def time_series_rmsd(u, arguments, verbose=False) -> List[float]:
    rmsds = []
    bb = u.select_atoms("backbone")

    J = len(u.trajectory)

    for t, traj in enumerate(u.trajectory):
        if t == 0:
            REF = bb.positions.copy()

            if arguments.ReferenceMean:
                ref_coordinates = u.trajectory.timeseries(asel=bb).mean(axis=1)
            else:
                ref_coordinates = u.trajectory.timeseries(asel=bb)[:, 0, :]

            # Make a reference structure (need to reshape into a
            # 1-frame "trajectory").
            ref = mda.Merge(bb).load_new(
                ref_coordinates[:, None, :],
                order="afc"
            )

            align.AlignTraj(
                u,
                ref,
                select="protein and name CA",
                in_memory=True
            ).run()

            # need to write the trajectory to
            # disk for PMDA 0.3.0 (see issue #15)
            with mda.Writer("rmsfit.xtc", n_atoms=u.atoms.n_atoms) as W:
                for ts in u.trajectory:
                    W.write(u.atoms)
        else:
            if verbose:
                print(f"{t} of {J}")

            v = rms.rmsd(REF, bb.positions)
            rmsds.append(v)

    return rmsds


def time_series_rmsf(u, end_pct=100) -> List[float]:
    sel = "protein and name CA"
    c_alphas = u.select_atoms(sel)
    return rms.RMSF(c_alphas).run().rmsf


def plotq(matrix):
    plt.imshow(matrix.dist_matrix, cmap='viridis')
    plt.xlabel('Frame')
    plt.ylabel('Frame')
    plt.colorbar(label=r'RMSD ($\AA$)')

    plt.show()


def main():
    arguments = parse_arguments()
    analyzeMD(arguments)


if __name__ == "__main__":
    main()
