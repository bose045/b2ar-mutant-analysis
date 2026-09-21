#!/usr/bin/env python3
"""Associate regional interhelical distance shifts with the response score S.

Each TM is represented by the centroid of three C-alpha atoms in its
extracellular, middle, or intracellular third/window, matching the Figure 4
analysis. All 21 TM pairs are measured in every retained frame in datasets A
and B. Ag-minus-IAg shifts are correlated with S across the seven mutants.
"""

from __future__ import annotations

from example_label import annotate_example
import math
import re
from itertools import combinations
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.io import netcdf_file
from scipy.signal import find_peaks
from scipy.stats import spearmanr


HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"

SYSTEMS = ["WT", "P60F", "P183H", "P228A", "P60M", "P183A", "P183F", "P228M"]
MUTANTS = SYSTEMS[1:]
TM_REGIONS = {
    "TM1": (2, 36), "TM2": (39, 68), "TM3": (77, 108), "TM4": (119, 142),
    "TM5": (169, 197), "TM6": (209, 238), "TM7": (245, 268),
}
TM_NAMES = list(TM_REGIONS)
PAIRS = list(combinations(TM_NAMES, 2))
PLANES = ["extracellular", "middle", "intracellular"]
ASCENDING = {"TM1", "TM3", "TM5", "TM7"}
N_BINS = 60


def basic_window(start, end, position, k=3):
    if position == "middle":
        midpoint = start + (end - start) // 2
        left = max(start, min(midpoint - k // 2, end - k + 1))
        return tuple(range(left, left + k))
    if position == "first":
        return tuple(range(start, start + k))
    return tuple(range(end - k + 1, end + 1))


def residues_for(tm, plane):
    start, end = TM_REGIONS[tm]
    if plane == "middle":
        return basic_window(start, end, "middle")
    if tm in ASCENDING:
        numeric = "first" if plane == "extracellular" else "last"
    else:
        numeric = "last" if plane == "extracellular" else "first"
    return basic_window(start, end, numeric)


WINDOWS = {plane: {tm: residues_for(tm, plane) for tm in TM_NAMES} for plane in PLANES}














def dominant_peak(values, edges):
    heights, _ = np.histogram(values, bins=edges, density=True)
    peaks, _ = find_peaks(heights, prominence=0.01)
    index = peaks[np.argmax(heights[peaks])] if len(peaks) else int(np.argmax(heights))
    return float((edges[index] + edges[index + 1]) / 2)


def summaries(raw, dataset):
    rows = []
    for system in SYSTEMS:
        for plane in PLANES:
            for tm1, tm2 in PAIRS:
                key = f"{plane}_{tm1}_{tm2}"
                active = raw[f"{system}_Active_{key}"]; inactive = raw[f"{system}_Inactive_{key}"]
                lo=float(min(active.min(),inactive.min()));hi=float(max(active.max(),inactive.max()))
                edges = np.linspace(lo if hi>lo else lo-.5,hi if hi>lo else hi+.5,N_BINS+1)
                rows.append({
                    "dataset": dataset, "system": system, "region": plane, "helix_pair": f"{tm1}-{tm2}",
                    "active_mean_A": active.mean(), "inactive_mean_A": inactive.mean(),
                    "mean_shift_A": active.mean() - inactive.mean(),
                    "median_shift_A": np.median(active) - np.median(inactive),
                    "peak_shift_A": dominant_peak(active, edges) - dominant_peak(inactive, edges),
                    "active_frames": len(active), "inactive_frames": len(inactive),
                })
    return pd.DataFrame(rows)


def associations(summary):
    rows = []
    for dataset in ("A", "B"):
        for plane in PLANES:
            for tm1, tm2 in PAIRS:
                pair = f"{tm1}-{tm2}"
                sub = summary[(summary.dataset == dataset) & (summary.region == plane) &
                              (summary.helix_pair == pair)].set_index("system").loc[MUTANTS]
                row = {"dataset": dataset, "region": plane, "helix_pair": pair}
                score = np.array([SCORES[dataset][m] for m in MUTANTS])
                for metric in ("mean_shift_A", "median_shift_A", "peak_shift_A"):
                    values=sub[metric].to_numpy()
                    rho,p = spearmanr(score,values) if np.ptp(values)>0 and np.ptp(score)>0 else (np.nan,np.nan)
                    row[f"{metric}_rho"] = rho; row[f"{metric}_p"] = p
                rows.append(row)
    return pd.DataFrame(rows)








def paired_matrix_plot_with_sign_disagreement_x(assoc):
    """Show dataset-specific median rho, replacing cross-dataset sign mismatches by x."""
    fig, axes = plt.subplots(2, 3, figsize=(15.5, 9.2))
    panel_labels = (("(a)", "(b)", "(c)"), ("(d)", "(e)", "(f)"))
    for col, region in enumerate(PLANES):
        region_data = assoc[assoc.region == region]
        lookup_a = region_data[region_data.dataset == "A"].set_index("helix_pair")
        lookup_b = region_data[region_data.dataset == "B"].set_index("helix_pair")
        disagreements = {}
        for pair in lookup_a.index:
            rho_a = float(lookup_a.loc[pair, "median_shift_A_rho"])
            rho_b = float(lookup_b.loc[pair, "median_shift_A_rho"])
            # A zero correlation has no direction and therefore does not satisfy
            # the same-sign requirement when the other dataset is nonzero.
            disagreements[pair] = np.isfinite(rho_a) and np.isfinite(rho_b) and np.sign(rho_a) != np.sign(rho_b)

        for row, dataset in enumerate(("A", "B")):
            ax = axes[row, col]
            matrix = np.full((7, 7), np.nan)
            crossed = []
            source = lookup_a if dataset == "A" else lookup_b
            for pair, rec in source.iterrows():
                tm1, tm2 = pair.split("-")
                i, j = TM_NAMES.index(tm1), TM_NAMES.index(tm2)
                if disagreements[pair]:
                    crossed.append((j, i))
                else:
                    matrix[j, i] = float(rec["median_shift_A_rho"])
            im = ax.imshow(np.ma.masked_invalid(matrix), cmap="RdBu_r", vmin=-1, vmax=1)
            ax.set_xticks(range(7), TM_NAMES, rotation=45, ha="right", fontsize=10)
            ax.set_yticks(range(7), TM_NAMES, fontsize=10)
            if row == 0:
                ax.set_title(region.capitalize(), fontsize=14, fontweight="bold", pad=10)
            ax.text(-.10, 1.04, panel_labels[row][col], transform=ax.transAxes,
                    ha="left", va="bottom", fontsize=14, fontweight="bold",
                    clip_on=False)
            for i in range(7):
                for j in range(7):
                    if np.isfinite(matrix[i, j]):
                        value = matrix[i, j]
                        ax.text(j, i, f"{value:+.2f}", ha="center", va="center", fontsize=9,
                                color="white" if abs(value) > .55 else "black")
            for i, j in crossed:
                ax.text(j, i, "×", ha="center", va="center", fontsize=15,
                        color="#666666", fontweight="bold")

    axes[0, 0].set_ylabel("Dataset A", fontsize=14, fontweight="bold", labelpad=12)
    axes[1, 0].set_ylabel("Dataset B", fontsize=14, fontweight="bold", labelpad=12)
    cbar = fig.colorbar(im, ax=axes, fraction=.022, pad=.025)
    cbar.set_label("Spearman ρ with S", fontsize=12)
    cbar.ax.tick_params(labelsize=10)
    fig.suptitle("Regional interhelical Ag−IAg median-distance shifts associated with S",
                 fontsize=17, fontweight="bold", y=.975)
    fig.text(.5, .025,
             "× indicates that the median-based correlations do not have the same sign in datasets A and B.",
             ha="center", fontsize=11)
    fig.subplots_adjust(left=.08, right=.90, top=.88, bottom=.12, wspace=.22, hspace=.28)
    annotate_example(fig);fig.savefig(RESULTS / "regional_interhelix_S_association_median_sign_disagreement_x.png",
                dpi=600, bbox_inches="tight", facecolor="white")
    annotate_example(fig);fig.savefig(RESULTS / "regional_interhelix_S_association_median_sign_disagreement_x.svg",
                bbox_inches="tight", facecolor="white")
    annotate_example(fig);fig.savefig(RESULTS / "regional_interhelix_S_association_median_sign_disagreement_x.pdf",
                bbox_inches="tight", facecolor="white")
    plt.close(fig)





# Supplied by input_root/regional_scores.json at runtime.
SCORES = {}
