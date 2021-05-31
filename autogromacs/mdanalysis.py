
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
    parser.add_argument("-d", dest='AutoDetect')

    parser.add_argument("-t", dest='TrajSuffix', default="")
    parser.add_argument("-M", dest='DoMatrix', action="store_true")
    parser.add_argument("-T", dest='DoTimeseries', action="store_true")
    parser.add_argument("-w", dest='WriteOutput', action="store_true")

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


def index_label(l: int) -> str:
    common = {
        0: "initial",
        -1: "last"
    }
    try:
        return common[l]
    except KeyError:
        return str(l)


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
    NB = int(re.findall(r"\d+", name)[0])
    if NB == 0:
        return "Original"

    return f"Variação #{NB}"


def RMSDStudy(us, unames):
    POS = []

    traj_idx = [0, -1]
    ATOM_ID = "name CA"
    ATOM_ID = "backbone"

    Data = []
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
    if arguments.AutoDetect:
        SimulationPrefixes = autodetect_files(arguments.AutoDetect)

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

    print("Select all or some? (input comma separated numbers and dash separated intervals)")

    q = input(">")

    q = [v.strip() for v in q.split(",")]

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
            OutputPrefixes.append(SimulationPrefixes[prefix_idx -1])

    for prefix in OutputPrefixes:
        print('\t' + prefix)

    return OutputPrefixes


def build_filepath(base: str, specifiers: List[str], arguments) -> Optional[str]:
    if arguments.WriteOutput:
        if arguments.AutoDetect:
            base = arguments.AutoDetect
            if arguments.FilePrefix:
                base += "+"
        else:
            base = "analysis"
        return base + "_" + "_".join(specifiers) + ".png"


def load_universe(SimulationPrefix, arguments):

    U = MDAnalysis.Universe(
        SimulationPrefix + ".gro",
        SimulationPrefix + arguments.TrajSuffix + ".trr"
    )

    aligner = align.AlignTraj(
        U,
        U,
        select='name CA',
        in_memory=True
    ).run()

    bb = U.select_atoms('protein and backbone')

    return U


def analyzeMD(arguments):

    SimulationPrefixes = loadSimulationPrefixes(arguments)
    SimulationPrefixes = selectSimulationPrefixes(SimulationPrefixes)
    #us = map(lambda sp: load_universe(sp, arguments), SimulationPrefixes)
    # can access via segid (4AKE) and atom name
    # we take the first atom named N and the last atom named C
    #nterm = u.select_atoms('segid 4AKE and name N')[0]
    #cterm = u.select_atoms('segid 4AKE and name C')[-1]

    base_filepath = arguments.WriteOutput if arguments.WriteOutput else None

    print("Data loading done.")

    if False:
        matrix = diffusionmap.DistanceMatrix(us[0], select='name CA').run()

    if arguments.DoMatrix:
        print("DEPRECATED.")
        sys.exit(1)

        print("Processing pairwise RMSD matrix.")
        POS = RMSDStudy(us, SimulationPrefixes)
        labels = [w.label for w in POS]
        RMSDS = pairwise_rmsds(POS)

        matrix_filepath = None
        if base_filepath is not None:
            matrix_filepath = base_filepath + "_matrix.png"

        show_matrix(RMSDS, labels, matrix_filepath)

    if arguments.DoTimeseries:
        print("Processing timeseries RMSD plots.")
        labels = []
        rmsd_series = []
        rmsf_series = []

        for SP in SimulationPrefixes:
            u = load_universe(SP, arguments)
            labels.append(get_label(u))
            rmsd_series.append(time_series_rmsd(u))
            rmsf_series.append(time_series_rmsf(u))

            del u

        show_rms_series(
            rmsd_series,
            labels,
            build_filepath(base_filepath, ["tsp", "rmsd"], arguments),
            "RMSD"
        )

        show_rms_series_monolithic(
            rmsd_series, labels,
            build_filepath(base_filepath, ["tsmono", "rmsd"], arguments),
            "RMSD"
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

    if False:
        for i, ts in enumerate(us[0].trajectory):
            # iterate through all frames
            #r = cterm.position - nterm.position
            # end-to-end vector from atom positions
            print(dir(ts))
            print(ts.positions)
            print(i)
    #r = bb.position
    #d = numpy.linalg.norm(r)  # end-to-end distance
    #rgyr = bb.radius_of_gyration()  # method of AtomGroup
    #print("frame = {0}: d = {1} A, Rgyr = {2} A".format(
    #      ts.frame, d, rgyr))


def show_matrix(results, labels, filepath: Union[str, None]):
    fig, ax = plt.subplots()

    im = ax.imshow(results, cmap='viridis')

    fig.colorbar(im, ax=ax, label=r'RMSD ($\AA$)')
    # We want to show all ticks...
    U = np.arange(len(labels))
    ax.set_yticks(U)
    plt.xticks(range(len(results)), labels, rotation='vertical')
    #ax.set_xticks(U, rotation='vertical')
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

    #ax.set_title(mode)
    ax.set_xlabel(XL)
    ax.set_ylabel(YL)

    ax.legend(labels)
    plt.tight_layout()

    if filepath is not None:
        plt.savefig(filepath)
    else:
        plt.show()


def show_rms_series(
        rms_series: List[List[float]],
        labels: List[str],
        filepath: Union[str, None],
        mode: str):

    X = len(labels)
    ncols = 3
    nrows = round(np.ceil(X / ncols))

    assert(ncols * nrows >= X)
    fig, ax = plt.subplots(nrows, ncols)

    axk = ax.ravel()

    for i, (vals, label) in enumerate(zip(rms_series, labels)):
        Xa = vals
        axk[i].plot(range(len(Xa)), Xa, "b-")

        axk[i].set_title(label)

        YL = r"Distância ($\AA$)"
        if mode == "RMSD":
            XL = "Frame"
        elif mode == "RMSF":
            XL = "Residue"
        else:
            exit(1)

        axk[i].set_xlabel(XL)
        axk[i].set_ylabel(YL)

    # plt.title(mode)
    plt.tight_layout()

    if filepath is not None:
        plt.savefig(filepath)
    else:
        plt.show()


def time_series_rmsd(u, verbose=False) -> List[float]:
    rmsds = []
    rmsfs = []
    bb = u.select_atoms("backbone")

    J = len(u.trajectory)

    for t, traj in enumerate(u.trajectory):
        if t == 0:
            REF = bb.positions.copy()

            #FREF = rms.RMSF(bb).run().rmsf
            ref_coordinates = u.trajectory.timeseries(asel=bb).mean(axis=1)

            # Make a reference structure (need to reshape into a
            # 1-frame "trajectory").
            ref = mda.Merge(bb).load_new(
                ref_coordinates[:, None, :],
                order="afc"
            )

            aligner = align.AlignTraj(
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


def time_series_rmsf(u) -> List[float]:
    bb = u.select_atoms("backbone")
    w = rms.RMSF(bb).run().rmsf
    print("RMSF")
    print(w.shape)
    print(w)

    return w


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

