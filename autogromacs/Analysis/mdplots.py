from typing import List, Union, Optional

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.units as munits

import numpy as np
import pandas as pd


def seaborn_theme():
    sns.set_theme(
        style="darkgrid",
        rc={
            'axes.facecolor': '#F2F2F2',
            'figure.facecolor': 'white'
        }
    )

    plt.rcParams["axes.labelsize"] = 15
    plt.rcParams["figure.dpi"] = 700


def _(_1, x, _2):
    return x


class ModeParameters():
    x_label = ""
    y_label = r"$\Delta$ RMSD ($\AA$)"
    enforce_ticks = True
    make_x = _

    def __init__(self, mode: str):

        if mode == "RMSDt":
            self.x_label = "Tempo (ns)"
            self.make_x = self.to_time_x
        elif mode == "RMSDf":
            self.x_label = "Frame"
        elif mode == "RMSF":
            self.x_label = "Resíduo"
            self.y_label = r"$\Delta$ RMSF ($\AA$)"
        elif mode == "PCA":
            self.x_label = "Componente"
            self.y_label = "Variância Acumulada"
        elif mode == "SASA":
            self.x_label = "Frame"
            self.y_label = r"SASA ($\AA$²)"
            self.enforce_ticks = False
        else:
            raise Exception("Unknown plot identifier.")

    def to_time_x(self, X, t):
        return frames_to_time(X, t)


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

    execute_output_plot(filepath)


def show_rms_series_monolithic(
        rms_series: List[List[float]],
        labels: List[str],
        total_times: List[float],
        filepath: Union[str, None],
        mode: str):

    fig, ax = plt.subplots()

    fig.set_figwidth(9.6)

    mode_parameters = ModeParameters(mode)

    values = np.array(rms_series)
    if values.ndim > 2:
        print(f"Incompatible data to plot for {filepath}.")
        return

    for i, Xa in enumerate(rms_series):
        x = mode_parameters.make_x(range(len(Xa)), total_times[i])
        ax.plot(x, Xa)

    # ax.set_title(mode)
    ax.set_xlabel(mode_parameters.x_label)
    ax.set_ylabel(mode_parameters.y_label)

    ax.legend(labels)

    plt.tight_layout()

    execute_output_plot(filepath)


def enforce_ax_ticks(ax, TICK_MAX, TICK_INTERVAL):

    YTICKS = list(range(0, 100, TICK_INTERVAL))

    # if Y_MAX < YTICK_INTERVAL:
    #     YTICKS.append(Y_MAX - 0.1)

    ax.set_yticks(sorted(YTICKS))

    ax.set_ylim(bottom=0, top=max(TICK_MAX, 1.05 * TICK_INTERVAL))


def hide_ax_ticks(ax):
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.set_facecolor('#ffffff')

    ax.tick_params(
        labelcolor='w',
        top=False,
        bottom=False,
        left=False,
        right=False
    )


def show_rms_series_stacked(
        rms_series: Union[List[List[float]], List[List[List[float]]]],
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
    fig.set_figheight(2.5 + V * 0.20 * N)
    fig.set_figwidth(V)

    ax = fig.add_subplot(111)
    axv = fig.subplots(nrows, ncols)

    Values = np.array(rms_series)

    try:
        Y_MAX = np.max(Values)
    except ValueError as e:
        raise e

    try:
        axk = axv.ravel()
    except AttributeError:
        axk = [axv]

    N = Values.shape[-1]

    def plot_bar(ax, X, Y, style, color):
        ax.bar(X, Y, color=color)

    def plot_line(ax, X, Y, style, color):
        ax.plot(X, Y, style, color=color)

    stacked_fn = [plot_bar, plot_line]
    for i, (Y, label) in enumerate(zip(rms_series, labels)):

        mode_parameters = ModeParameters(mode)

        X: Union[List[float], List[int]] = mode_parameters.make_x(
            list(range(N)),
            total_times[i]
        )
        print(Values.ndim)
        if Values.ndim == 3:
            colors = ["black", "orange"]
            styles = ["-", "--"]
            INITIALIZED = False

            for sY, plot_fn, style, color in zip(Y, stacked_fn, styles, colors):
                if INITIALIZED:
                    cax = axk[i].twinx()

                    hide_ax_ticks(cax)
                else:
                    cax = axk[i]
                try:
                    plot_fn(cax, X, sY, style, color)
                except munits.ConversionError as e:
                    print(X)
                    print(sY)
                    assert len(X) == len(sY)
                    raise(e)

                if mode_parameters.enforce_ticks:
                    enforce_ax_ticks(cax, round(max(sY)), round(max(sY)))
                INITIALIZED = True

        elif Values.ndim == 2:
            axk[i].plot(X, Y, "-", color="black")
            if mode_parameters.enforce_ticks:
                enforce_ax_ticks(axk[i], Y_MAX, 5)
        else:
            print(f"{Values.ndim}")
            raise Exception(f"Unexpected values shape of {Values.shape}")

        axk[i].set_ylabel(label, fontsize=12)
        axk[i].yaxis.set_label_position("right")

        axk[i].grid(b=False, axis='x')
        axk[i].tick_params(bottom=False)

        if i + 1 < len(labels):
            axk[i].set(xlabel=None)

    # fig.text(0.5, 0.01, XL, ha='center')
    # fig.text(0.00, 0.5, YL, va='center', rotation='vertical')

    hide_ax_ticks(ax)
    ax.set_xlabel(mode_parameters.x_label)
    ax.set_ylabel(mode_parameters.y_label)

    plt.subplots_adjust(hspace=0.05)

    execute_output_plot(filepath)


def execute_output_plot(filepath: Optional[str]) -> None:
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


def data_sanity(rms_values):
    pass


def plot_ramachandran(rama_analysis, label: str, output_file: str):
    fig, ax = plt.subplots(figsize=plt.figaspect(1))
    rama_analysis.plot(ax=ax, color='k', marker='o', s=5, ref=True)
    #fig.title(label)
    fig.tight_layout()
    fig.savefig(output_file)


seaborn_theme()
