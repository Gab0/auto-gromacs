 
from typing import List, Union

import seaborn as sns
import matplotlib.pyplot as plt

import numpy as np

sns.set_theme(style="darkgrid")
plt.rcParams["axes.labelsize"] = 15


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

    def to_time_x(X):
        return frames_to_time(X, 64)

    def _(x):
        return x

    YL = r"Distância ($\AA$)"
    if mode == "RMSDt":
        XL = "Tempo (ns)"
        make_x = to_time_x
    elif mode == "RMSDf":
        XL = "Frame"
        make_x = _

    elif mode == "RMSF":
        XL = "Residue"
    elif mode == "PCA":
        XL = "Component"
        YL = "% Variance"
    else:
        raise Exception("Unknown plot identifier.")

    for i, Xa in enumerate(rms_series):
        ax.plot(make_x(range(len(Xa))), Xa)

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
    ncols = 1
    nrows = round(np.ceil(N / ncols))

    assert ncols * nrows >= N

    fig = plt.figure()
    ax = fig.add_subplot(111)
    axv = fig.subplots(nrows, ncols)

    try:
        axk = axv.ravel()
    except AttributeError:
        axk = [axv]

    for i, (vals, label) in enumerate(zip(rms_series, labels)):

        Y = vals
        X: Union[List[float], List[int]] = list(range(len(Y)))

        YL = r"Distância ($\AA$)"
        if mode == "RMSDf":
            XL = "Frame"
        elif mode == "RMSDt":
            XL = "Tempo (ns)"
            X = frames_to_time(X, 64)
        elif mode == "RMSF":
            XL = "Residue"
        else:
            exit(1)

        axk[i].plot(X, Y, "b-")
        axk[i].set_title(label)

    # fig.text(0.5, 0.01, XL, ha='center')
    # fig.text(0.00, 0.5, YL, va='center', rotation='vertical')
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.set_facecolor('#ffffff')

    ax.tick_params(labelcolor='w', top=False,
                   bottom=False, left=False, right=False)

    ax.set_xlabel(XL)
    ax.set_ylabel(YL)
    # plt.title(mode)
    plt.tight_layout()

    if filepath is not None:
        plt.savefig(filepath)
    else:
        plt.show()
