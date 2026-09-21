#!/usr/bin/env python
"""Plot the submitted Figure 5 analysis using Dataset-B contact frequencies."""

from __future__ import annotations

from example_label import annotate_example
import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import networkx as nx
import numpy as np
import pandas as pd

from network_plotting_common import (
    CONTACT_MARKERS,
    KEYS,
    TM_ORDER,
    draw_network,
    legend_handles,
    residue_label,
    winning_forms,
)


def strongest_residue_edges(edges: pd.DataFrame) -> pd.DataFrame:
    """Collapse parallel contact types to the largest-|F| edge for cluster display."""
    ordered = edges.sort_values("abs_functional_shift_msgm", ascending=False)
    return ordered.drop_duplicates(["res1_num", "res2_num"]).copy()


def build_residue_graph(edges: pd.DataFrame) -> nx.Graph:
    graph = nx.Graph()
    for _, row in strongest_residue_edges(edges).iterrows():
        first = residue_label(row["res1_name"], row["res1_num"])
        second = residue_label(row["res2_name"], row["res2_num"])
        graph.add_edge(first, second, **row.to_dict())
    return graph


def draw_clusters(
    axis: plt.Axes,
    edges: pd.DataFrame,
    panel_label: str = "(b)",
    font_scale: float = 1.0,
    node_scale: float = 1.0,
    panel_font_scale: float | None = None,
    layout_scale: float = 1.0,
) -> pd.DataFrame:
    panel_font_scale = font_scale if panel_font_scale is None else panel_font_scale
    graph = build_residue_graph(edges)
    components = [sorted(component) for component in nx.connected_components(graph) if len(component) >= 4]
    components.sort(key=lambda component: (-len(component), component))

    records = []
    if not components:
        axis.text(0.5, 0.5, "No connected components with ≥4 residues",
                  transform=axis.transAxes, ha="center", va="center")
        axis.axis("off")
        return pd.DataFrame()

    widths = [max(1.0, np.sqrt(len(component))) for component in components]
    centers = np.cumsum([0.0] + widths[:-1]) + np.asarray(widths) / 2
    centers = 0.08 + 0.84 * centers / (centers[-1] + widths[-1] / 2)
    global_min = float(edges.abs_functional_shift_msgm.min())
    global_max = float(edges.abs_functional_shift_msgm.max())

    for component_index, (component, center, width) in enumerate(zip(components, centers, widths), start=1):
        subgraph = graph.subgraph(component).copy()
        local = nx.spring_layout(subgraph, seed=42 + component_index, iterations=400)
        xs = np.asarray([local[node][0] for node in component])
        ys = np.asarray([local[node][1] for node in component])
        xscale = 0.34 * layout_scale * width / max(sum(widths), 1.0)
        yscale = 0.34 * layout_scale
        positions = {
            node: (center + xscale * local[node][0], 0.50 + yscale * local[node][1])
            for node in component
        }
        for first, second, data in subgraph.edges(data=True):
            score = float(data["functional_shift_msgm"])
            magnitude = float(data["abs_functional_shift_msgm"])
            edge_width = 1.2 + 4.8 * (magnitude - global_min) / max(global_max - global_min, 1e-12)
            color = "#e31a1c" if score > 0 else "#314cff"
            style = "solid" if winning_forms(pd.Series(data)) else (0, (5, 3))
            x = [positions[first][0], positions[second][0]]
            y = [positions[first][1], positions[second][1]]
            axis.plot(x, y, color=color, linewidth=edge_width, linestyle=style,
                      alpha=0.78, zorder=1)
            records.append({
                "component": component_index,
                "component_size": len(component),
                "residue_1": first,
                "residue_2": second,
                "contact_type": data["contact_type"],
                "functional_shift_msgm": score,
            })
        coords = np.asarray([positions[node] for node in component])
        axis.scatter(coords[:, 0], coords[:, 1], s=300 * node_scale, facecolor="#dddddd",
                     edgecolor="#555555", linewidth=0.7, zorder=3)
        for node in component:
            axis.text(*positions[node], node, fontsize=6.2 * font_scale,
                      ha="center", va="center", zorder=4)

    axis.text(-0.02, 1.02, panel_label, transform=axis.transAxes,
              fontsize=17 * panel_font_scale, fontweight="bold", va="top")
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)
    axis.axis("off")
    return pd.DataFrame(records)
