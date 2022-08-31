from glob import glob
from os.path import basename
from typing import Sequence

import intervaltree
import matplotlib.pyplot as plt
import numpy as np

# plt.rcParams["backend"] = "svg"

# from scipy.signal import savgol_filter

# TODO add "PIL" to the conda env


def max_min_score(scores_dir: str) -> tuple[float, float]:
    """
    for the BWA score get the max/min score to use for scaling
    """
    max_score: float = 0
    min_score: float = 0
    beds = glob(f"{scores_dir}/*.bed")
    for bed in beds:
        with open(bed) as b:
            for line in b:
                v = float(line.rstrip().split()[3])
                if v > max_score:
                    max_score = v
                if v < min_score:
                    min_score = v
    return max_score, min_score


def min_max_scale(x: float, min_x: float, max_x: float) -> float:
    """
    Scale x between 0 and 1
    """
    return (x - min_x) / (max_x - min_x)


def smooth_data(
    x: Sequence,
    window_size: int,
) -> np.ndarray:
    """
    Smooth the data using a moving average filter
    """
    return np.convolve(x, np.ones((window_size,)) / window_size, mode="same")


def plot_scores(
    scores_dir: str,
    scores_type: str,
    line_width: float = 0.5,
    xlim: tuple[float, float] = (0.0, 30000.0),
    ylim: tuple[float, float] = (0.0, 1.0),
    smoothing_window: int = 10,
    figsize: tuple[int, int] = (10, 10),
) -> None:
    """
    Plot the scores accross the virus genome of the various query sequences.
    The name of each query sequence is the name of the file containing the scores.

    Inputs:
    - scores_dir: path to the directory containing the beds with scores
    - scores_type: {"bwa" | "similiarity"}
    """

    # TODO let max/min score be optional(float|None)
    if scores_type == "bwa":
        max_score, min_score = max_min_score(scores_dir)
    else:
        max_score = -1
        min_score = -1

    beds = glob(f"{scores_dir}/*.bed")
    plt.figure(figsize=figsize)
    labels: list[str] = []
    for bed in sorted(beds):
        interval_bins = intervaltree.IntervalTree()
        for i in range(0, 30_000, 50):
            interval_bins.addi(i, i + 50, [])

        bed_intervals = intervaltree.IntervalTree()
        with open(bed) as b:
            for line in b:
                A = line.rstrip().split()
                s, e = map(int, A[1:3])
                if max_score > 0:
                    v = min_max_scale(float(A[3]), min_score, max_score)
                else:
                    v = float(A[3])
                bed_intervals.addi(s, e, v)

        for i in bed_intervals:
            ovlps = interval_bins.overlap(i)
            for o in ovlps:
                o.data.append(i.data)

        bins = sorted(
            [(i.begin, i.end, np.mean(i.data)) for i in interval_bins if i.data],
            key=lambda x: x[0],
        )

        labels.append(basename(bed).split("_")[0])
        if scores_type == "similarity":
            plt.plot(
                [x[0] for x in bins],
                # savgol_filter([1 - np.sqrt(x[2]) / 2 for x in bins], 51, 3),
                smooth_data([1 - np.sqrt(x[2]) / 2 for x in bins], smoothing_window),
                label=labels[-1],
                linewidth=line_width,
            )
        elif scores_type == "bwa":
            plt.plot(
                [x[0] for x in bins],
                # savgol_filter([x[2] for x in bins], 51, 3),
                smooth_data([x[2] for x in bins], smoothing_window),
                label=labels[-1],
                linewidth=line_width,
            )
        else:
            raise ValueError('scores_type must be "similarity" or "bwa"')

    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.xlabel("genomic position")
    plt.ylabel("similarity score")
    plt.gca().spines['bottom'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    leg = plt.legend(loc="center left", bbox_to_anchor=(1, 0.5), labels=sorted(labels))
    # plt.legend(
    #     loc="upper center",
    #     bbox_to_anchor=(0.5, 1.5),
    #     ncol=3,
    #     fancybox=True,
    #     shadow=True,
    # )
    # plt.legend()
    for legobj in leg.legendHandles:
        legobj.set_linewidth(2.0)
    plt.show()
