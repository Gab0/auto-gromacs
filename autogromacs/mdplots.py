from typing import List, Union, Optional

import seaborn as sns
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd

sns.set_theme(style="darkgrid")
plt.rcParams["axes.labelsize"] = 15
plt.rcParams["figure.dpi"] = 700


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

    execute_output(filepath)


def show_rms_series_monolithic(
        rms_series: List[List[float]],
        labels: List[str],
        total_times: List[float],
        filepath: Union[str, None],
        mode: str):

    fig, ax = plt.subplots()

    fig.set_figwidth(9.6)

    def to_time_x(X, t):
        return frames_to_time(X, t)

    def _(X, t):
        return X

    YL = r"Distância ($\AA$)"
    make_x = _
    if mode == "RMSDt":
        XL = "Tempo (ns)"
        make_x = to_time_x
    elif mode == "RMSDf":
        XL = "Frame"
    elif mode == "RMSF":
        XL = "Resíduo"
    elif mode == "PCA":
        XL = "Componente"
        YL = "Variância Acumulada"
    else:
        raise Exception("Unknown plot identifier.")

    for i, Xa in enumerate(rms_series):
        ax.plot(make_x(range(len(Xa)), total_times[i]), Xa)

    # ax.set_title(mode)
    ax.set_xlabel(XL)
    ax.set_ylabel(YL)

    ax.legend(labels)
    plt.tight_layout()

    execute_output(filepath)


def execute_output(filepath: Optional[str]) -> None:
    if filepath is not None:
        plt.savefig(filepath)
    else:
        plt.show()


def frames_to_time(frames: Union[List[float], List[int]],
                   total_time: float) -> List[float]:
    total_frames = frames[-1]

    def to_t(v: Union[int, float]) -> float:
        return v / total_frames * total_time

    return list(map(to_t, frames))


def write_series(
        series_data: List[List[float]],
        labels: List[str],
        filepath: str):

    data = {}
    for d, l in zip(series_data, labels):
        data[l] = d
    df = pd.DataFrame(data, columns=labels)
    df.to_csv(filepath)


def show_rms_series(
        rms_series: List[List[float]],
        labels: List[str],
        total_times: List[float],
        filepath: Union[str, None],
        mode: str):

    N = len(labels)
    ncols = 1
    nrows = round(np.ceil(N / ncols))

    assert ncols * nrows >= N

    fig = plt.figure()

    V = 6
    fig.set_figheight(1 + V * 0.2 * N)
    fig.set_figwidth(V)

    ax = fig.add_subplot(111)
    axv = fig.subplots(nrows, ncols)

    Values = np.array(rms_series)
    Y_MAX = np.max(Values)
    Y_MIN = np.min(Values)

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
            X = frames_to_time(X, total_times[i])
        elif mode == "RMSF":
            XL = "Residue"
        else:
            exit(1)

        axk[i].plot(X, Y, "b-")

        axk[i].set_ylabel(label, fontsize=12)
        axk[i].yaxis.set_label_position("right")

        axk[i].set_ylim(bottom=Y_MIN, top=Y_MAX)
        if i + 1 < len(labels):
            axk[i].set(xlabel=None)

        axk[i].grid(b=False, axis='x')
        axk[i].tick_params(bottom=False)

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

    plt.subplots_adjust(hspace=0.05)

    if filepath is not None:
        plt.savefig(filepath)
    else:
        plt.show()
