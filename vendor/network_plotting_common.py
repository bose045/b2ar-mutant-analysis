#!/usr/bin/env python
"""Generate lenient and cross-criterion robust functional-shift networks for SI."""

from __future__ import annotations

from example_label import annotate_example
import argparse
from collections import defaultdict
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from matplotlib.path import Path as MplPath
from matplotlib.patches import PathPatch
import numpy as np
import pandas as pd


TM_REGIONS = {
    "TM1": (1, 36), "TM2": (39, 68), "TM3": (77, 108), "TM4": (119, 142),
    "TM5": (169, 197), "TM6": (209, 238), "TM7": (245, 268), "H8": (270, 280),
}
TM_ORDER = list(TM_REGIONS)
TM_COLORS = {
    "TM1": "#1f77b4", "TM2": "#ff7f0e", "TM3": "#d62728", "TM4": "#8c564b",
    "TM5": "#e377c2", "TM6": "#bcbd22", "TM7": "#17becf", "H8": "#7f7f7f",
}
CONTACT_MARKERS = {
    "hbbb": ("H-bond BB", "o"), "hbsb": ("H-bond SB", "^"),
    "hbss": ("H-bond SS", "s"), "hp": ("Hydrophobic", "p"),
    "vdw": ("van der Waals", "X"), "wb": ("Water bridge", "D"),
    "pc": ("Cation–π", "*"), "ps": ("π-stack", "P"),
    "ts": ("T-stack", "h"), "sb": ("Salt bridge", "d"),
}
KEYS = ["res1_num", "res2_num", "contact_type"]


def assign_tm(residue: int) -> str | None:
    for name, (start, end) in TM_REGIONS.items():
        if start <= int(residue) <= end:
            return name
    return None


def residue_label(name: str, number: int) -> str:
    return f"{name}-{int(number)}"


def residue_number(label: str) -> int:
    return int(label.rsplit("-", 1)[1])


def make_robust_table(analysis: Path) -> pd.DataFrame:
    criteria = ("strict", "default", "lenient")
    tables = {
        criterion: pd.read_csv(analysis / f"display_edges_F0.3_{criterion}.csv")
        for criterion in criteria
    }
    base_columns = KEYS + ["res1_name", "res2_name"]
    merged = tables["default"][base_columns].copy()
    for criterion in criteria:
        columns = KEYS + [
            "functional_shift_msgm", "msgm_act", "msgm_inact",
            "maj_sign_act", "maj_sign_inact", "n_agree_act", "n_agree_inact",
        ]
        renamed = tables[criterion][columns].rename(
            columns={column: f"{column}_{criterion}" for column in columns if column not in KEYS}
        )
        merged = merged.merge(renamed, on=KEYS, how="inner")

    score_columns = [f"functional_shift_msgm_{criterion}" for criterion in criteria]
    score_signs = np.sign(merged[score_columns].to_numpy(float))
    sign_consistent = np.all(score_signs == score_signs[:, [0]], axis=1)
    if not np.all(sign_consistent):
        raise ValueError("At least one intersection edge changes functional-shift sign across criteria")

    robust = merged.loc[sign_consistent].copy()
    robust["functional_shift_msgm"] = robust[score_columns].median(axis=1)
    robust["abs_functional_shift_msgm"] = robust["functional_shift_msgm"].abs()
    for column in ("msgm_act", "msgm_inact"):
        robust[column] = robust[[f"{column}_{criterion}" for criterion in criteria]].median(axis=1)
    robust["maj_sign_act"] = np.sign(robust["msgm_act"])
    robust["maj_sign_inact"] = np.sign(robust["msgm_inact"])
    robust["n_agree_act"] = robust[[f"n_agree_act_{c}" for c in criteria]].min(axis=1)
    robust["n_agree_inact"] = robust[[f"n_agree_inact_{c}" for c in criteria]].min(axis=1)
    robust["criterion"] = "robust_all_three"
    robust["robust_definition"] = (
        "Passes |F|>=0.3 and all eligibility rules under strict, default, and lenient criteria; "
        "plotted score is the median across criteria"
    )
    return robust.sort_values(KEYS).reset_index(drop=True)


def winning_forms(row: pd.Series) -> bool:
    winning = row["msgm_act"] if row["functional_shift_msgm"] > 0 else row["msgm_inact"]
    return float(winning) > 0


def sign_flip(row: pd.Series) -> bool:
    first, second = float(row["maj_sign_act"]), float(row["maj_sign_inact"])
    return first != 0 and second != 0 and first != second


def layout_nodes(edges: pd.DataFrame) -> tuple[list[str], dict[str, tuple[float, float]], dict[str, float]]:
    tm_nodes = {tm: set() for tm in TM_ORDER}
    for row in edges.itertuples(index=False):
        for name, number in ((row.res1_name, row.res1_num), (row.res2_name, row.res2_num)):
            tm = assign_tm(number)
            if tm:
                tm_nodes[tm].add(residue_label(name, number))
    ordered = {tm: sorted(nodes, key=residue_number) for tm, nodes in tm_nodes.items()}
    present = [tm for tm in TM_ORDER if ordered[tm]]
    nodes = [node for tm in present for node in ordered[tm]]
    radius, gap = 10.0, np.deg2rad(12.0)
    step = (2 * np.pi - gap * max(0, len(present) - 1)) / len(nodes)
    position: dict[str, tuple[float, float]] = {}
    angle: dict[str, float] = {}
    current = np.deg2rad(90.0)
    for tm in present:
        for node in ordered[tm]:
            angle[node] = current
            position[node] = (radius * np.cos(current), radius * np.sin(current))
            current -= step
        if tm != present[-1]:
            current -= gap
    return nodes, position, angle


def draw_network(
    axis: plt.Axes,
    edges: pd.DataFrame,
    title: str,
    width_min: float,
    width_max: float,
    panel_label: str | None = None,
    label_scores: bool = False,
    font_scale: float = 1.0,
    residue_font_scale: float | None = None,
    score_font_scale: float | None = None,
    tm_font_scale: float | None = None,
    panel_font_scale: float | None = None,
    residue_label_radius: float = 10.8,
    tm_label_radius: float = 12.6,
    axis_limit: float = 14.2,
) -> tuple[list[str], list[str]]:
    residue_font_scale = font_scale if residue_font_scale is None else residue_font_scale
    score_font_scale = font_scale if score_font_scale is None else score_font_scale
    tm_font_scale = font_scale if tm_font_scale is None else tm_font_scale
    panel_font_scale = font_scale if panel_font_scale is None else panel_font_scale
    nodes, position, angle = layout_nodes(edges)
    present_tms = [tm for tm in TM_ORDER if any(assign_tm(residue_number(n)) == tm for n in nodes)]
    radius = 10.0

    def edge_width(value: float) -> float:
        return 0.9 + 4.1 * (value - width_min) / max(width_max - width_min, 1e-12)

    grouped: dict[tuple[str, str], list[pd.Series]] = defaultdict(list)
    for _, row in edges.iterrows():
        first = residue_label(row["res1_name"], row["res1_num"])
        second = residue_label(row["res2_name"], row["res2_num"])
        grouped[tuple(sorted((first, second)))].append(row)

    for (first, second), rows in grouped.items():
        p0, p2 = np.asarray(position[first]), np.asarray(position[second])
        chord = p2 - p0
        length = float(np.hypot(*chord)) or 1.0
        perpendicular = np.array([-chord[1] / length, chord[0] / length])
        midpoint = 0.5 * (p0 + p2)
        offsets = (np.arange(len(rows)) - (len(rows) - 1) / 2) * 0.55
        for offset, row in zip(offsets, rows):
            color = "#d62728" if row["functional_shift_msgm"] > 0 else "#1f4ed8"
            line_style = "solid" if winning_forms(row) else (0, (5, 3))
            line_width = edge_width(float(row["abs_functional_shift_msgm"]))
            control = 0.52 * midpoint + offset * perpendicular
            path = MplPath([tuple(p0), tuple(control), tuple(p2)],
                           [MplPath.MOVETO, MplPath.CURVE3, MplPath.CURVE3])
            axis.add_patch(PathPatch(path, facecolor="none", edgecolor="#555555",
                                     linewidth=line_width + 1.2, alpha=0.32, zorder=1))
            axis.add_patch(PathPatch(path, facecolor="none", edgecolor=color,
                                     linewidth=line_width, linestyle=line_style,
                                     alpha=0.78, zorder=2))
            marker = CONTACT_MARKERS.get(row["contact_type"], ("Other", "x"))[1]
            point = 0.25 * p0 + 0.5 * control + 0.25 * p2
            axis.scatter(*point, marker=marker, s=34, facecolor=color,
                         edgecolor="black", linewidth=0.45, zorder=4)
            if label_scores:
                label_point = point - 0.28 * perpendicular
                axis.text(*label_point, f"{float(row['functional_shift_msgm']):.2f}",
                          color=color, fontsize=6.2 * score_font_scale, fontweight="bold",
                          ha="center", va="center", zorder=8)
            if sign_flip(row):
                axis.scatter(*(point + 0.10 * perpendicular), marker="*", s=30,
                             facecolor="#ffdf00", edgecolor="black", linewidth=0.45, zorder=5)

    coordinates = np.array([position[node] for node in nodes])
    node_colors = [TM_COLORS[assign_tm(residue_number(node))] for node in nodes]
    axis.scatter(coordinates[:, 0], coordinates[:, 1], s=175, c=node_colors,
                 edgecolors="black", linewidths=0.65, zorder=6)

    for node in nodes:
        theta = angle[node]
        x, y = residue_label_radius * np.cos(theta), residue_label_radius * np.sin(theta)
        rotation, alignment = np.degrees(theta), "left"
        if np.cos(theta) < 0:
            rotation += 180
            alignment = "right"
        axis.text(x, y, node, fontsize=7.6 * residue_font_scale, rotation=rotation,
                  rotation_mode="anchor", ha=alignment, va="center", zorder=7)

    for tm in present_tms:
        tm_angles = [angle[node] for node in nodes if assign_tm(residue_number(node)) == tm]
        middle = float(np.mean(tm_angles))
        axis.text(tm_label_radius * np.cos(middle), tm_label_radius * np.sin(middle), tm,
                  fontsize=10.5 * tm_font_scale, fontweight="bold", color=TM_COLORS[tm],
                  ha="center", va="center", zorder=7)

    axis.set_title(title, fontsize=13 * font_scale, pad=16)
    if panel_label:
        axis.text(-0.03, 1.03, panel_label, transform=axis.transAxes,
                  fontsize=14 * panel_font_scale, fontweight="bold", va="top")
    axis.set_aspect("equal")
    axis.set_xlim(-axis_limit, axis_limit)
    axis.set_ylim(-axis_limit, axis_limit)
    axis.axis("off")
    return nodes, present_tms


def legend_handles(
    edges: pd.DataFrame,
    present_tms: list[str],
    marker_scale: float = 1.0,
) -> list[Line2D | mpatches.Patch]:
    handles: list[Line2D | mpatches.Patch] = [
        Line2D([], [], linestyle="none", label="TM region"),
        *[mpatches.Patch(color=TM_COLORS[tm], label=tm) for tm in present_tms],
        Line2D([], [], linestyle="none", label="Contact type"),
    ]
    for contact_type in CONTACT_MARKERS:
        if contact_type in set(edges["contact_type"]):
            label, marker = CONTACT_MARKERS[contact_type]
            handles.append(Line2D([], [], linestyle="none", marker=marker,
                                  markerfacecolor="#777777", markeredgecolor="black",
                                  markersize=6 * marker_scale, label=label))
    handles.extend([
        Line2D([], [], linestyle="none", label="Functional-shift direction"),
        Line2D([], [], color="#d62728", linewidth=2, label="Activating indicator (+F)"),
        Line2D([], [], color="#1f4ed8", linewidth=2, label="Inactivating indicator (−F)"),
        Line2D([], [], linestyle="none", label="Contact change"),
        Line2D([], [], color="#666666", linewidth=2, linestyle="solid", label="Forms in active state"),
        Line2D([], [], color="#666666", linewidth=2, linestyle=(0, (5, 3)), label="Breaks in active state"),
        Line2D([], [], linestyle="none", marker="*", markerfacecolor="#ffdf00",
               markeredgecolor="black", markersize=7 * marker_scale,
               label="Opposite signs between classes"),
    ])
    return handles
