"""
plot_SLC16A3_expression.py
---------------
Reads SLC16A3_results_breast_cell_lines.xls (Results sheet) and produces
a publication-quality bar chart of relative SLC16A3 expression (2^-ddCt)
normalised to geometric mean of ACTB and RNA18SN5, with MCF10A as reference.

Usage:
    python plot_SLC16A3_expression.py

Output:
    SLC16A3_expression.png  (300 dpi, for paper)
    SLC16A3_expression.pdf  (vector, for submission)

Requirements:
    pip install pandas xlrd matplotlib numpy scipy
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats
import warnings
warnings.filterwarnings("ignore")

# ── 1. Load data ──────────────────────────────────────────────────────────────
FILE = "SLC16A3_results_breast_cell_lines.xls"

df_raw = pd.read_excel(FILE, sheet_name="Results", engine="xlrd", header=None)

# Row index 7 is the true header
df = df_raw.iloc[7:].copy()
df.columns = df_raw.iloc[7].values
df = df.iloc[1:].reset_index(drop=True)

# Keep only plate well rows
df = df[df["Well"].astype(str).str.match(r"^[A-Z]\d+$")].copy()
df["Cт"] = pd.to_numeric(df["Cт"], errors="coerce")
df = df[["Well", "Sample Name", "Target Name", "Cт"]].dropna()

# Rename targets and samples to standard names
target_map = {"SLC": "SLC16A3", "18S": "RNA18SN5", "BACT": "ACTB"}
sample_map = {
    "MCF10A KONTROLA": "MCF10A",
    "MCF7":            "MCF7",
    "MDA231":          "MDA-MB-231",
    "MDA468":          "MDA-MB-468",
}
df["Target Name"] = df["Target Name"].map(target_map)
df["Sample Name"] = df["Sample Name"].map(sample_map)
df = df.dropna(subset=["Sample Name", "Target Name"])

# ── 2. Replicates: assign replicate index (1,2,3) within each sample+target ──
df = df.sort_values(["Sample Name", "Target Name", "Well"]).reset_index(drop=True)
df["rep"] = df.groupby(["Sample Name", "Target Name"]).cumcount() + 1

# ── 3. Pivot so each row = one replicate of one sample ───────────────────────
pivot = df.pivot_table(index=["Sample Name", "rep"],
                       columns="Target Name",
                       values="Cт",
                       aggfunc="first").reset_index()
pivot = pivot.dropna(subset=["SLC16A3", "RNA18SN5", "ACTB"])

# ── 4. ΔCt = SLC16A3 − geometric mean of housekeepers ────────────────────────
# Geometric mean of Ct values on log scale = arithmetic mean
pivot["geo_ref"] = (pivot["RNA18SN5"] + pivot["ACTB"]) / 2
pivot["dCt"]     = pivot["SLC16A3"] - pivot["geo_ref"]

# ── 5. ΔΔCt relative to MCF10A mean ΔCt ────────────────────────────────────
ref_mean = pivot[pivot["Sample Name"] == "MCF10A"]["dCt"].mean()
pivot["ddCt"]     = pivot["dCt"] - ref_mean
pivot["rel_expr"] = 2 ** (-pivot["ddCt"])

# ── 6. Summary stats ──────────────────────────────────────────────────────────
order = ["MCF10A", "MCF7", "MDA-MB-231", "MDA-MB-468"]

summary = (
    pivot.groupby("Sample Name")["rel_expr"]
    .agg(
        mean="mean",
        sem=lambda x: x.std(ddof=1) / np.sqrt(len(x)),
    )
    .reindex(order)
)

means  = summary["mean"].values.astype(float)
sems   = summary["sem"].values.astype(float)

# Individual replicate values for dot overlay
rep_vals = {c: pivot[pivot["Sample Name"] == c]["rel_expr"].values for c in order}

# ── 7. Statistics vs MCF10A (one-sample t-test, µ = 1) ───────────────────────
def sig_label(p):
    if p < 0.001: return "***"
    elif p < 0.01:  return "**"
    elif p < 0.05:  return "*"
    else:           return "ns"

pvals = {}
for cell in order[1:]:
    _, p = stats.ttest_1samp(rep_vals[cell], popmean=1.0)
    pvals[cell] = p

# ── 8. Plot ───────────────────────────────────────────────────────────────────
COLORS = {
    "MCF10A":    "#4A7FB5",
    "MCF7":      "#C0392B",
    "MDA-MB-231":"#E67E22",
    "MDA-MB-468":"#27AE60",
}
EDGE = "#1a1a1a"
GREY = "#7f7f7f"

fig, ax = plt.subplots(figsize=(6.5, 5.2))
fig.patch.set_facecolor("white")
ax.set_facecolor("white")

x = np.arange(len(order))
bar_w = 0.52

# Bars
ax.bar(x, means, width=bar_w,
       color=[COLORS[c] for c in order],
       edgecolor=EDGE, linewidth=0.9,
       zorder=3)

# Error bars separately so they sit on top
ax.errorbar(x, means, yerr=sems,
            fmt="none", elinewidth=1.3, capsize=5, capthick=1.3,
            ecolor=EDGE, zorder=5)

# Individual data points
rng = np.random.default_rng(0)
for i, cell in enumerate(order):
    pts = rep_vals[cell]
    jitter = rng.uniform(-0.07, 0.07, size=len(pts))
    ax.scatter(x[i] + jitter, pts,
               color="white", edgecolors=EDGE,
               linewidths=0.9, s=32, zorder=6)

# Reference line at y = 1
ax.axhline(1, color=GREY, linewidth=0.85, linestyle="--", zorder=2)

# Significance brackets
def draw_bracket(ax, x1, x2, y, label, dy):
    ax.plot([x1, x1, x2, x2],
            [y, y + dy, y + dy, y],
            color=EDGE, lw=0.9)
    ax.text((x1 + x2) / 2, y + dy + 0.012,
            label, ha="center", va="bottom", fontsize=9.5)

top = float(np.nanmax(means + sems)) * 1.15
step = top * 0.22
dy   = top * 0.07

for i, cell in enumerate(order[1:], start=1):
    draw_bracket(ax, 0, i, top + step * (i - 1), sig_label(pvals[cell]), dy)

ymax = top + step * (len(order) - 1) + dy * 2.5
ax.set_ylim(-0.04, ymax)

# Axes cosmetics
ax.set_xticks(x)
ax.set_xticklabels(order, fontsize=10.5, fontstyle="italic")
ax.set_ylabel(
    "Relative SLC16A3 expression\n(2$^{-\\Delta\\Delta Ct}$, relative to MCF10A)",
    fontsize=10
)
ax.set_title("SLC16A3 mRNA expression in breast cell lines",
             fontsize=11.5, fontweight="bold", pad=10)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_linewidth(0.9)
ax.spines["bottom"].set_linewidth(0.9)
ax.tick_params(axis="both", length=3, width=0.9, labelsize=9)
ax.yaxis.set_minor_locator(plt.MultipleLocator(0.1))

plt.tight_layout()
plt.savefig("SLC16A3_expression.png", dpi=300, bbox_inches="tight")
plt.savefig("SLC16A3_expression.pdf", bbox_inches="tight")
print("Saved: SLC16A3_expression.png and SLC16A3_expression.pdf")
