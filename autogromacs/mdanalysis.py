
from typing import List
import argparse
import sys
import numpy as np
import warnings
import MDAnalysis
from MDAnalysis.analysis import diffusionmap, align, rms

import numpy.linalg

import matplotlib.pyplot as plt


warnings.filterwarnings("ignore", category=DeprecationWarning)

def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument("-m", dest='FilePrefix', nargs="*")


    return parser.parse_args()


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


def clean_universe_prefix(p):
    try:
        p = p.split("/")[-2]
    except IndexError:
        pass

    p = p.replace("work", "")
    p = p.strip("_")

    return p


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


def analyzeMD(arguments):
    print("File prefixes used:")
    for prefix in arguments.FilePrefix:
        print("\t" + prefix)

    us = [
        MDAnalysis.Universe(FP + ".gro", FP + ".trr")
        for FP in arguments.FilePrefix
    ]

    # can access via segid (4AKE) and atom name
    # we take the first atom named N and the last atom named C
    #nterm = u.select_atoms('segid 4AKE and name N')[0]
    #cterm = u.select_atoms('segid 4AKE and name C')[-1]

    for U in us:
        aligner = align.AlignTraj(U,
                                  U,
                                  select='name CA',
                                  in_memory=True).run()
    bb = [
        u.select_atoms('protein and backbone')
        for u in us
    ] # a selection (AtomGroup)


    matrix = diffusionmap.DistanceMatrix(us[0], select='name CA').run()

    if True:
        print("Processing pairwise RMSD matrix.")
        POS = RMSDStudy(us, arguments.FilePrefix)
        labels = [w.label for w in POS]
        RMSDS = pairwise_rmsds(POS)
        print(RMSDS)
        show_matrix(RMSDS, labels)

    if True:
        print("Processing timeseries RMSD plots.")
        labels = [
            clean_universe_prefix(f)
            for f in arguments.FilePrefix
        ]
        series = list(map(time_series_rms, us))
        show_rmsd_series(series, labels)

    sys.exit(0)
    for i, ts in enumerate(us[0].trajectory):     # iterate through all frames
        #r = cterm.position - nterm.position # end-to-end vector from atom positions
        print(dir(ts))
        print(ts.positions)
        print(i)
    #r = bb.position
    #d = numpy.linalg.norm(r)  # end-to-end distance
    #rgyr = bb.radius_of_gyration()  # method of AtomGroup
    #print("frame = {0}: d = {1} A, Rgyr = {2} A".format(
    #      ts.frame, d, rgyr))


def show_matrix(results, labels):
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

    plt.show()


def show_rmsd_series(rmsd_series: List[List[List[float]]],
                     labels: List[str]):
    X = len(labels)
    ncols = 3
    nrows = round(np.ceil(X / ncols))

    assert(ncols * nrows >= X)
    fig, ax = plt.subplots(nrows, ncols)

    axk = ax.ravel()

    for i, (vals, label) in enumerate(zip(rmsd_series, labels)):
        Xa = vals[0]
        Xb = vals[1]

        axk[i].plot(range(len(Xa)), Xa, "b-")
        axk[i].plot(range(len(Xb)), Xb, "r-")
        axk[i].set_title(label)

    plt.tight_layout()
    plt.show()


def time_series_rms(u, verbose=True):
    rmsds = []
    rmsfs = []
    bb = u.select_atoms("backbone")
    J = len(u.trajectory)

    for t, traj in enumerate(u.trajectory):
        if t == 0:
            REF = bb.positions.copy()
            #FREF = rms.RMSF(bb).run().rmsf
        else:
            if verbose:
                print(f"{t} of {J}")

            v = rms.rmsd(REF, bb.positions)
            rmsds.append(v)

            #w = np.mean(rms.RMSF(bb).run().rmsf - FREF)
            #rmsfs.append(w)

    return [ rmsds, [] ]


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

