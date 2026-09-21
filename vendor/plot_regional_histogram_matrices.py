#!/usr/bin/env python3
"""Plot matched regional interhelical-distance histogram matrices for datasets A/B.

Uses the cached framewise centroid distances produced by
``analyze_interhelix_regions.py``; trajectory data are not reread.
"""

from __future__ import annotations

from example_label import annotate_example
import math
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch
import numpy as np


HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
OUT = RESULTS / "regional_histogram_matrices"
HIGHLIGHT_OUT = OUT / "opposite_consensus_highlighted"
SYSTEMS = ["WT", "P60F", "P183H", "P228A", "P60M", "P183A", "P183F", "P228M"]
REGIONS = ["extracellular", "middle", "intracellular"]
TMS = [f"TM{i}" for i in range(1, 8)]
ACTIVE_COLOR = "#0072B2"
INACTIVE_COLOR = "#D55E00"
ACTIVE_MEDIAN_COLOR = "#003F63"
INACTIVE_MEDIAN_COLOR = "#7A3000"
AL_ARROW_COLOR = "#159B35"
IL_ARROW_COLOR = "#9A159B"
OPPOSITE_HALO_COLOR = "#FFD400"
SHIFT_THRESHOLD = 0.25
S_THRESHOLD = 1.0


def fd_edges(active: np.ndarray, inactive: np.ndarray) -> np.ndarray:
    combined = np.concatenate((active, inactive))
    q25, q75 = np.percentile(combined, [25, 75])
    width = 2.0 * (q75 - q25) / np.cbrt(len(combined))
    if not np.isfinite(width) or width <= 0 or combined.max() == combined.min():
        return np.linspace(combined.min() - 0.5, combined.max() + 0.5, 41)
    count = max(15, min(100, int(math.ceil((combined.max() - combined.min()) / width))))
    return np.linspace(combined.min(), combined.max(), count + 1)




def functional_class(dataset: str, system: str):
    if system == "WT":
        return None
    score = SCORES[dataset][system]
    if score >= S_THRESHOLD:
        return "AL"
    if score <= -S_THRESHOLD:
        return "IL"
    return None


def pair_key(system: str, state: str, region: str, tm1: str, tm2: str):
    return f"{system}_{state}_{region}_{tm1}_{tm2}"


def median_shift(cache, system: str, region: str, tm1: str, tm2: str):
    active = cache[pair_key(system, "Active", region, tm1, tm2)]
    inactive = cache[pair_key(system, "Inactive", region, tm1, tm2)]
    return float(np.median(active) - np.median(inactive))


def build_class_consensus(caches):
    """Return >=3-of-5 direction consensus for pooled S-classified observations."""
    consensus = {}
    records = []
    for region in REGIONS:
        for i in range(1, 7):
            for j in range(i):
                tm1, tm2 = TMS[j], TMS[i]
                for klass in ("AL", "IL"):
                    observations = []
                    for dataset in ("A", "B"):
                        for system in SYSTEMS[1:]:
                            if functional_class(dataset, system) == klass:
                                shift = median_shift(caches[dataset], system, region, tm1, tm2)
                                observations.append((dataset, system, shift))
                    positive = sum(shift >= SHIFT_THRESHOLD for _, _, shift in observations)
                    negative = sum(shift <= -SHIFT_THRESHOLD for _, _, shift in observations)
                    if positive>=3 and negative>=3:
                        raise ValueError('Both signs satisfy fixed three-observation consensus; review changed class sizes')
                    sign = 1 if positive >= 3 else (-1 if negative >= 3 else 0)
                    consensus[(region, tm1, tm2, klass)] = sign
                    records.append({
                        "region": region, "helix_pair": f"{tm1}-{tm2}", "class": klass,
                        "n_observations": len(observations), "positive_resolved": positive,
                        "negative_resolved": negative, "consensus_sign": sign,
                        "consensus_direction": "Ag separation" if sign > 0 else
                                               ("Ag contraction" if sign < 0 else "none"),
                    })
    return consensus, records


def consensus_relation(consensus, region: str, tm1: str, tm2: str, klass):
    """Return ``opposite``, ``exclusive``, or ``none`` for a class consensus."""
    if klass not in ("AL", "IL"):
        return "none"
    other = "IL" if klass == "AL" else "AL"
    sign = consensus.get((region, tm1, tm2, klass), 0)
    other_sign = consensus.get((region, tm1, tm2, other), 0)
    if sign and other_sign == -sign:
        return "opposite"
    if sign and other_sign == 0:
        return "exclusive"
    return "none"


def informative_consensus_sign(consensus, region: str, tm1: str, tm2: str, klass):
    """Return class sign for opposite-class or class-exclusive consensuses.

    Same-direction consensuses shared by AL and IL are intentionally suppressed.
    """
    sign = consensus.get((region, tm1, tm2, klass), 0)
    return sign if consensus_relation(consensus, region, tm1, tm2, klass) != "none" else 0


def highlight_opposite(arrow, consensus, region, tm1, tm2, klass):
    """Add a gold halo when AL and IL have opposite consensus signs."""
    if consensus_relation(consensus, region, tm1, tm2, klass) == "opposite":
        arrow.set_path_effects([
            pe.Stroke(linewidth=5.8, foreground=OPPOSITE_HALO_COLOR),
            pe.Normal(),
        ])
    return arrow






def make_region_comparison_figure(caches, region: str, consensus):
    """Create one figure comparing datasets A and B for a single region."""
    fig = plt.figure(figsize=(18, 18.2), facecolor="white")
    outer = fig.add_gridspec(2, 1, left=.035, right=.985, bottom=.035, top=.915,
                             hspace=.12)

    for dataset_index, dataset in enumerate(("A", "B")):
        cache = caches[dataset]
        section = outer[dataset_index].subgridspec(
            3, 4, height_ratios=(.11, 1, 1), wspace=.20, hspace=.17
        )
        header = fig.add_subplot(section[0, :])
        header.axis("off")
        header.text(.002, .5, f"({'ab'[dataset_index]})", ha="left", va="center",
                    fontsize=19, fontweight="bold")
        header.text(.5, .5, f"Dataset {dataset}", ha="center", va="center",
                    fontsize=20, fontweight="bold")

        for panel_index, system in enumerate(SYSTEMS):
            slot = section[panel_index // 4 + 1, panel_index % 4]
            grid = slot.subgridspec(7, 7, wspace=.045, hspace=.045)
            for row in range(7):
                for col in range(7):
                    ax = fig.add_subplot(grid[row, col])
                    if row < col:
                        ax.axis("off")
                        continue
                    if row == col:
                        ax.axis("off")
                        ax.text(.5, .5, TMS[row], ha="center", va="center",
                                fontsize=9, fontweight="bold")
                        continue

                    tm1, tm2 = TMS[col], TMS[row]
                    active = cache[pair_key(system, "Active", region, tm1, tm2)]
                    inactive = cache[pair_key(system, "Inactive", region, tm1, tm2)]
                    edges = fd_edges(active, inactive)
                    ax.hist(inactive, bins=edges, density=True, color=INACTIVE_COLOR,
                            alpha=.88, linewidth=0)
                    ax.hist(active, bins=edges, density=True, color=ACTIVE_COLOR,
                            alpha=.88, linewidth=0)
                    med_active = float(np.median(active))
                    med_inactive = float(np.median(inactive))
                    ax.axvline(med_active, color=ACTIVE_MEDIAN_COLOR, lw=1.3, zorder=8)
                    ax.axvline(med_inactive, color=INACTIVE_MEDIAN_COLOR, lw=1.3, zorder=8)

                    klass = functional_class(dataset, system)
                    shift = med_active - med_inactive
                    shift_sign = (1 if shift >= SHIFT_THRESHOLD else
                                  (-1 if shift <= -SHIFT_THRESHOLD else 0))
                    class_sign = informative_consensus_sign(
                        consensus, region, tm1, tm2, klass
                    )
                    if shift_sign and shift_sign == class_sign:
                        ymax = ax.get_ylim()[1]
                        arrow = FancyArrowPatch(
                            (med_inactive, .83 * ymax), (med_active, .83 * ymax),
                            arrowstyle="-|>", mutation_scale=15, linewidth=2.35,
                            color=AL_ARROW_COLOR if klass == "AL" else IL_ARROW_COLOR,
                            zorder=12,
                        )
                        highlight_opposite(arrow, consensus, region, tm1, tm2, klass)
                        ax.add_patch(arrow)
                    ax.set_yticks([])
                    ax.tick_params(axis="x", labelsize=6.8, length=1.8, pad=1)
                    if row != 6:
                        ax.set_xticklabels([])
                    else:
                        ticks = np.linspace(edges[0], edges[-1], 3)[1:]
                        ax.set_xticks(ticks)
                        ax.set_xticklabels([f"{x:.0f}" for x in ticks])
                    for spine in ax.spines.values():
                        spine.set_color("#C9C9C9")
                        spine.set_linewidth(.55)

            box = slot.get_position(fig)
            fig.text(box.x0 + .55 * box.width, box.y1 - .002, system,
                     ha="center", va="top", fontsize=17, fontweight="bold")

    handles = [
        Line2D([0], [0], color=ACTIVE_COLOR, lw=8, label="Ag-bound"),
        Line2D([0], [0], color=INACTIVE_COLOR, lw=8, label="IAg-bound"),
        Line2D([0], [0], color=ACTIVE_MEDIAN_COLOR, lw=1.8, label="Ag median"),
        Line2D([0], [0], color=INACTIVE_MEDIAN_COLOR, lw=1.8, label="IAg median"),
        Line2D([0], [0], color=AL_ARROW_COLOR, lw=2.8, marker=">",
               markersize=8, label="Activation-like direction"),
        Line2D([0], [0], color=IL_ARROW_COLOR, lw=2.8, marker=">",
               markersize=8, label="Inactivation-like direction"),
        Line2D([0], [0], color=OPPOSITE_HALO_COLOR, lw=5,
               label="Gold halo: opposite class consensuses"),
    ]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(.5, .953),
               frameon=False, fontsize=11.5, handlelength=2.6, ncol=7,
               columnspacing=1.3)
    fig.suptitle(f"{region.capitalize()} interhelical-distance distributions",
                 fontsize=23, fontweight="bold", y=.986)
    fig.text(.5, .009, "Distance between three-Cα regional centroids (Å)",
             ha="center", fontsize=13)
    return fig





# Supplied by input_root/regional_scores.json at runtime.
SCORES = {}
