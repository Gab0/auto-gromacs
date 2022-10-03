from typing import List, Union, Optional, Callable, Any

import itertools
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.units as munits

import numpy as np
import pandas as pd


def seaborn_theme():
    """All plots will use this seaborn theme."""
    sns.set_theme(
        style="darkgrid",
        rc={
            'axes.facecolor': '#F2F2F2',
            'figure.facecolor': 'white',
            'font.family': 'sans-serif',
            'font.sans-serif': ['Montserrat', 'Verdana']
        }
    )

    plt.rcParams["axes.labelsize"] = 15
    plt.rcParams["figure.dpi"] = 700


def _(_0, _1, _2):
    """Helper function, just returns the second argument."""
    return _1


class ModeParameters():
    """Decodes the mode string into various plot parameter sets."""
    label_time = "Tempo (ns)"
    x_label = label_time
    y_label = r"$\Delta$ RMSD ($\AA$)"
    y_label2 = None

    enforce_ticks = True
    make_x: Callable[[Any, List[float], Any], List[float]] = _
    moving_average = False

    barline = False

    def __init__(self, mode: Optional[str]):
        if mode is None:
            pass
        elif mode == "RMSDt":
            self.make_x = self.to_time_x
        elif mode == "RMSDf":
            self.x_label = "Frame"
        elif mode == "RMSF":
            self.x_label = "Resíduo"
            self.y_label = r"$\Delta$ RMSF ($\AA$)"
        elif mode == "PCA":
            self.x_label = "Componente"
            self.y_label = "Variância Acumulada"
            self.barline = True
        elif mode == "SASA":
            self.y_label = r"SASA ($\AA$²)"
            self.make_x = self.to_time_x
        elif mode == "RADGYR":
            self.x_label = "Frame"
            self.y_label = r"Raio de Giro ($\AA$)"
        elif mode == "NSECONDARY":
            self.y_label = r"$n$ de Resíduos em Estruturas Secundárias"
            self.make_x = self.to_time_x
            self.moving_average = True
        elif mode == "ANGLES":
            self.y_label = "Ângulo entre os domínios"
            self.make_x = self.to_time_x
        elif mode == "RMSDANGLES":
            self.y_label2 = "Ângulo entre os domínios"
            self.make_x = self.to_time_x
        elif mode == "RMSDRMSD":
            self.y_label2 = self.y_label
            self.make_x = self.to_time_x
        else:
            raise Exception(f"Unknown plot identifier: {mode}.")

    def to_time_x(self, x_values, _time):
        return frames_to_time(x_values, _time)


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


def enforce_ax_ticks(ax, tick_min, tick_max, tick_interval):

    YTICKS = list(range(int(tick_min),
                        int(tick_max),
                        int(tick_interval)))

    ax.set_yticks(sorted(YTICKS))

    margin_mod = 0.05 * (tick_max - tick_min)

    ax.set_ylim(
        bottom=tick_min - margin_mod,
        top=tick_max + margin_mod
    )


def hide_ax_ticks(ax):
    """Hide axis spines for all directions.."""
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


def determine_tick_interval(y_max, y_min) -> int:

    y_delta = int(y_max - y_min)
    optimal_interval = 5
    if y_delta < optimal_interval:
        return y_delta
    if y_max <= 25:
        return optimal_interval

    return int(y_delta // 5)


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
        Y_MIN = np.min(Values)
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

        plot_initialized = False
        # Plots with three dimensions such as PCA plots;
        if Values.ndim == 3 and mode_parameters.barline:
            colors = ["black", "orange"]
            styles = ["-", "--"]

            for sY, plot_fn, style, color in zip(Y, stacked_fn, styles, colors):
                if plot_initialized:
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
                    enforce_ax_ticks(cax, round(min(sY)), round(max(sY)), round(max(sY)))

                plot_initialized = True

        # Common 2D plots (lines);
        else:
            # Ensure correct shape for Y values,
            # considering the initial Y may have 2D or 3D data (with multiple lines per axis).
            # FIXME: Yeah this whole operation is confusing.
            check_Y = np.array(Y)
            if check_Y.ndim == 1:
                Y = check_Y.reshape(1, -1)

            colors = ["black", "dimgrey"]

            for aY, color in zip(Y, colors):
                print(aY)
                if plot_initialized:
                    cax = axk[i].twinx()

                    hide_ax_ticks(cax)
                else:
                    cax = axk[i]

                plot_initialized = True

                _X = X
                if len(X) != len(aY):
                    _X = adapt_length_x(X, len(aY))

                try:
                    cax.plot(_X, aY, "-", color=color)
                except ValueError:
                    print("aY ERROR!")
                    cax.plot(_X, aY[0], "-", color=color)

                if mode_parameters.enforce_ticks:
                    tick_interval = determine_tick_interval(Y_MAX, Y_MIN)
                    enforce_ax_ticks(cax, Y_MIN, Y_MAX, tick_interval)
                if mode_parameters.moving_average:
                    moving_average = pd.Series(aY).rolling(10).mean()
                    cax.plot(X, moving_average, "-", color="grey", linewidth=1.5)
            #else:
            #    print(f"Number of dimensions for value: {Values.ndim}")
            #    raise Exception(f"Unexpected values shape of {Values.shape}")

        if False:
            axk[i].set_ylabel(label, fontsize=12)
            axk[i].yaxis.set_label_position("right")
        else:
            axk[i].text(0.1, 0.75, label, fontsize=14, transform=axk[i].transAxes)

        axk[i].grid(b=False, axis='x')
        axk[i].tick_params(bottom=False)

        if i + 1 < len(labels):
            axk[i].set(xlabel=None)
            axk[i].get_xaxis().set_ticks([])
    # fig.text(0.5, 0.01, XL, ha='center')
    # fig.text(0.00, 0.5, YL, va='center', rotation='vertical')

    y_pad = calculate_pad(len(str(Y_MAX)))
    print(f"Y label pad: {y_pad}")
    hide_ax_ticks(ax)
    ax.set_xlabel(mode_parameters.x_label, labelpad=0.6)
    ax.set_ylabel(mode_parameters.y_label, labelpad=y_pad)

    if mode_parameters.y_label2 is not None:
        ax2 = ax.twinx()
        ax2.set_ylabel(mode_parameters.y_label2, labelpad=y_pad)

    plt.subplots_adjust(hspace=0.03)
    #plt.tight_layout()
    execute_output_plot(filepath)


def calculate_pad(n):
    return n * 4.0


def execute_output_plot(filepath: Optional[str]) -> None:
    if filepath is not None:
        plt.savefig(filepath, bbox_inches="tight")
    else:
        plt.show()


def frames_to_time(frames: Union[List[float], List[int]],
                   total_time: float) -> List[float]:
    total_frames = frames[-1]

    def to_t(v: Union[int, float]) -> float:
        return v / total_frames * total_time

    return list(map(to_t, frames))


def match_series_sampling(A, B):
    """
    Equalize the lengths of two time series,
    making them proportional along a line graph.
    """
    target_length = max([len(k) for k in [A, B]])

    def match_series(s):
        if len(s) == target_length:
            return s

    return (match_series(k) for k in (A, B))


def adapt_length_x(x, target_length):
    modifier = target_length / len(x)

    #return [v * modifier for v in x]
    return list(range(target_length))


def write_series(
        series_data: List[List[float]],
        labels: List[str],
        filepath: str):

    data = {}
    for d, l in zip(series_data, labels):
        data[l] = d
    df = pd.DataFrame(data, columns=labels)
    df.to_csv(filepath)


def check_data_sanity(rms_values) -> bool:
    pass


def plot_ramachandran(rama_analysis, label: str, output_file: str):
    fig, ax = plt.subplots(figsize=plt.figaspect(1))
    rama_analysis.plot(ax=ax, color='k', marker='o', s=5, ref=True)
    #fig.title(label)
    fig.tight_layout()
    fig.savefig(output_file)


seaborn_theme()
