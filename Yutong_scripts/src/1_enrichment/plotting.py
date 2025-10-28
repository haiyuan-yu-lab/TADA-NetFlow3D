import argparse
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

def _p_to_stars(p):
    if p is None or np.isnan(p):
        return ""
    if p < 0.0001: return "****"
    if p < 0.001: return "***"
    if p < 0.01: return "**"
    if p < 0.05: return "*"

def plot_logOR_horizontal(data, se_scale, annotate="stars",
                          add_values=False, group_gap=0.6, label_fontsize=12,
                          tick_fontsize=12, x_label="Log Odds Ratio",
                          title=None, fig_size=(7, 0), output=None):
    """
    Flexible plotting of horizontal log(OR) with error bars and p-value annotations.

    Parameters
    ----------
    data : dict
        Two accepted structures:
          A) Flat dict:
              { "EnhancerType": (logOR, SE, pval), ... }
          B) Nested dict (grouped):
              { "EnhancerType": {
                    "combined": (logOR, SE, pval),
                    "shared":   (logOR, SE, pval),
                    "unique":   (logOR, SE, pval)
                 }, ... }
        Values must be (logOR, SE, pval). SE is the standard error of log(OR).
    """
    # -----------------------
    # Normalize input to rows
    # -----------------------
    rows = []  # each: dict(group, name, logor, se, pval)
    is_grouped = any(isinstance(v, dict) for v in data.values())
    if is_grouped:
        for group, subdict in data.items():
            if not isinstance(subdict, dict):
                raise ValueError("Mixed data structures: some groups are not dicts.")
            items = list(subdict.items())
            for name, (logor, se, pval) in items:
                rows.append({"group": group, "name": name, "logor": logor, "se": se, "pval": pval})
    else:
        items = list(data.items())
        for group, tup in items:
            logor, se, pval = tup
            rows.append({"group": group, "name": None, "logor": logor, "se": se, "pval": pval})

    # -----------------------
    # Build y positions
    # -----------------------
    y_positions = []
    y_labels = []
    current_y = 0.0
    last_group = None

    if is_grouped:
        # Preserve original group order
        group_order = []
        seen = set()
        for r in rows:
            g = r["group"]
            if g not in seen:
                seen.add(g)
                group_order.append(g)
        # regroup rows by group to inject gaps and labels
        grouped = {g: [] for g in group_order}
        for r in rows:
            grouped[r["group"]].append(r)

        flat_rows = []
        for g in group_order:
            if flat_rows:
                current_y += group_gap   # add a gap between different top-level groups
            for r in grouped[g]:
                y_positions.append(current_y)
                # label: "Group (Sub)" if sub-name exists; otherwise just group
                if r["name"]:
                    y_labels.append(f"{r['group']} ({r['name']})")
                else:
                    y_labels.append(f"{r['group']}")
                flat_rows.append(r)
                current_y += 1.0
        rows = flat_rows
    else:
        # Simple flat case
        for r in rows:
            y_positions.append(current_y)
            y_labels.append(r["group"])
            current_y += 1.0

    y_positions = np.array(y_positions)

    # -----------------------
    # Collect x, xerr
    # -----------------------
    x = np.array([r["logor"] for r in rows], dtype=float)
    se = np.array([r["se"] for r in rows], dtype=float) * float(se_scale)
    pvals = [r["pval"] for r in rows]

    # -----------------------
    # Figure & axes
    # -----------------------
    n = len(rows)
    if fig_size[1] == 0:
        # heuristic: ~0.5 inch per row + top/bottom margins
        height = max(2.5, 0.55 * n + 1.3 + (0.15 if is_grouped else 0))
        fig_size = (fig_size[0], height)

    fig, ax = plt.subplots(figsize=fig_size)

    # ---- errorbar points ----
    _point_kwargs = dict(fmt='o', capsize=5, elinewidth=1.5)
    ax.errorbar(x, y_positions, xerr=se, **_point_kwargs)

    # ---- zero line ----
    _line_kwargs = dict(linestyle='--', color="grey", linewidth=1.0)
    ax.axvline(0.0, **_line_kwargs)

    # -----------------------
    # Annotations
    # -----------------------
    for xi, yi, s, p in zip(x, y_positions, se, pvals):
        if annotate == "stars":
            label = _p_to_stars(p)
        elif annotate == "pval":
            label = f"p={p:.2e}" if p is not None and not np.isnan(p) else ""
        else:
            label = ""
        if label:
            # offset to the right of the error bar
            ax.text(xi + s + 0.02, yi, label, va="center", fontsize=tick_fontsize)

        if add_values:
            ax.text(xi, yi, f"{xi:.2f}", va="bottom", ha="center", fontsize=tick_fontsize-1)

    # -----------------------
    # Cosmetics
    # -----------------------
    ax.set_yticks(y_positions)
    ax.set_yticklabels(y_labels, fontsize=label_fontsize)
    ax.set_ylim(y_positions.min() - 0.5, y_positions.max() + 0.5)
    ax.set_xlabel(x_label, fontsize=label_fontsize + 2)
    if title:
        ax.set_title(title, fontsize=label_fontsize + 2)

    ax.invert_yaxis()  # top row at top
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.tick_params(axis="x", labelsize=tick_fontsize)
    ax.tick_params(axis="y", labelsize=tick_fontsize)

    plt.tight_layout()
    if output:
        plt.savefig(output, bbox_inches="tight")
    plt.show()
    return fig, ax