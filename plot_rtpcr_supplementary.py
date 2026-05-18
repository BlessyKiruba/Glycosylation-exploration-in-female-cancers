"""
plot_rtpcr_supplementary.py
---------------------
Generates three supplementary figures for SLC16A3 RT-qPCR data:

  S1 — Raw Ct values for all targets across all cell lines
  S2 — ΔCt values per cell line (normalized, pre-reference)
  S3 — Intra-assay replicate consistency heatmap

Usage:
    python plot_rtpcr_supplementary.py

Output:
    FigS1_raw_ct.png / .pdf
    FigS2_delta_ct.png / .pdf
    FigS3_replicate_heatmap.png / .pdf

Requirements:
    pip install pandas xlrd matplotlib numpy seaborn
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import warnings
warnings.filterwarnings("ignore")

# ═══════════════════════════════════════════════════════════════════════════════
# 0. SHARED DATA LOADING & CALCULATION
# ═══════════════════════════════════════════════════════════════════════════════

FILE = "SLC16A3_results_breast_cell_lines.xls"

df_raw = pd.read_excel(FILE, sheet_name="Results", engine="xlrd", header=None)
df = df_raw.iloc[7:].copy()
df.columns = df_raw.iloc[7].values
df = df.iloc[1:].reset_index(drop=True)
df = df[df["Well"].astype(str).str.match(r"^[A-Z]\d+$")].copy()
df["Cт"] = pd.to_numeric(df["Cт"], errors="coerce")
df = df[["Well", "Sample Name", "Target Name", "Cт"]].dropna()

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
df = df.sort_values(["Sample Name", "Target Name", "Well"]).reset_index(drop=True)
df["rep"] = df.groupby(["Sample Name", "Target Name"]).cumcount() + 1

pivot = df.pivot_table(
    index=["Sample Name", "rep"],
    columns="Target Name",
    values="Cт",
    aggfunc="first"
).reset_index().dropna(subset=["SLC16A3", "RNA18SN5", "ACTB"])

pivot["geo_ref"] = (pivot["RNA18SN5"] + pivot["ACTB"]) / 2
pivot["dCt"]     = pivot["SLC16A3"] - pivot["geo_ref"]
ref_mean         = pivot[pivot["Sample Name"] == "MCF10A"]["dCt"].mean()
pivot["ddCt"]    = pivot["dCt"] - ref_mean
pivot["rel_expr"]= 2 ** (-pivot["ddCt"])

ORDER   = ["MCF10A", "MCF7", "MDA-MB-231", "MDA-MB-468"]
TARGETS = ["SLC16A3", "RNA18SN5", "ACTB"]
COLORS  = {
    "MCF10A":    "#4A7FB5",
    "MCF7":      "#C0392B",
    "MDA-MB-231":"#E67E22",
    "MDA-MB-468":"#27AE60",
}
TARGET_COLORS = {
    "SLC16A3":  "#2C3E50",
    "RNA18SN5": "#8E44AD",
    "ACTB":     "#16A085",
}
EDGE = "#1a1a1a"
GREY = "#7f7f7f"
rng  = np.random.default_rng(42)

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE S1 — Raw Ct values for all targets
# ═══════════════════════════════════════════════════════════════════════════════

fig1, axes = plt.subplots(1, 3, figsize=(12, 4.8), sharey=False)
fig1.patch.set_facecolor("white")

for ax, target in zip(axes, TARGETS):
    sub = df[df["Target Name"] == target]
    x   = np.arange(len(ORDER))

    for i, cell in enumerate(ORDER):
        vals = sub[sub["Sample Name"] == cell]["Cт"].values
        mean = vals.mean()
        sem  = vals.std(ddof=1) / np.sqrt(len(vals))

        ax.bar(i, mean, width=0.55,
               color=COLORS[cell], edgecolor=EDGE,
               linewidth=0.8, zorder=3, alpha=0.88)
        ax.errorbar(i, mean, yerr=sem,
                    fmt="none", elinewidth=1.2, capsize=4,
                    capthick=1.2, ecolor=EDGE, zorder=5)

        jitter = rng.uniform(-0.08, 0.08, size=len(vals))
        ax.scatter(i + jitter, vals,
                   color="white", edgecolors=EDGE,
                   linewidths=0.8, s=28, zorder=6)

    # Y axis range — zoom in to show variation clearly
    all_vals = sub["Cт"].dropna().values
    ypad = (all_vals.max() - all_vals.min()) * 0.55
    ax.set_ylim(all_vals.min() - ypad, all_vals.max() + ypad)

    ax.set_xticks(x)
    ax.set_xticklabels(ORDER, fontsize=9, fontstyle="italic")
    ax.set_title(f"$\\it{{{target}}}$", fontsize=11, fontweight="bold",
                 color=TARGET_COLORS[target])
    ax.set_ylabel("C$_t$ value" if target == "SLC16A3" else "", fontsize=9.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.8)
    ax.spines["bottom"].set_linewidth(0.8)
    ax.tick_params(length=3, width=0.8, labelsize=8.5)



plt.tight_layout(w_pad=3)
plt.savefig("FigS1_raw_ct.png", dpi=300, bbox_inches="tight")
plt.savefig("FigS1_raw_ct.pdf", bbox_inches="tight")
print("Saved FigS1_raw_ct.png / .pdf")
plt.close()

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE S2 — ΔCt values per cell line
# ═══════════════════════════════════════════════════════════════════════════════

fig2, ax2 = plt.subplots(figsize=(6, 4.8))
fig2.patch.set_facecolor("white")
ax2.set_facecolor("white")

x = np.arange(len(ORDER))

dct_summary = (
    pivot.groupby("Sample Name")["dCt"]
    .agg(mean="mean", sem=lambda v: v.std(ddof=1) / np.sqrt(len(v)))
    .reindex(ORDER)
)

for i, cell in enumerate(ORDER):
    mean = dct_summary.loc[cell, "mean"]
    sem  = dct_summary.loc[cell, "sem"]
    vals = pivot[pivot["Sample Name"] == cell]["dCt"].values

    ax2.bar(i, mean, width=0.55,
            color=COLORS[cell], edgecolor=EDGE,
            linewidth=0.8, zorder=3, alpha=0.88)
    ax2.errorbar(i, mean, yerr=sem,
                 fmt="none", elinewidth=1.2, capsize=4,
                 capthick=1.2, ecolor=EDGE, zorder=5)

    jitter = rng.uniform(-0.08, 0.08, size=len(vals))
    ax2.scatter(i + jitter, vals,
                color="white", edgecolors=EDGE,
                linewidths=0.8, s=32, zorder=6)

    # Annotate mean value
    ax2.text(i, mean + sem + 0.07, f"{mean:.2f}",
             ha="center", va="bottom", fontsize=8, color=EDGE)

ax2.set_xticks(x)
ax2.set_xticklabels(ORDER, fontsize=10.5, fontstyle="italic")
ax2.set_ylabel("ΔC$_t$  (SLC16A3 − geometric mean of ACTB, RNA18SN5)",
               fontsize=9.5)
ax2.set_title("Figure S2 — ΔC$_t$ values per cell line\n"
              "(higher ΔC$_t$ = lower relative SLC16A3 expression)",
              fontsize=11, fontweight="bold", pad=10)

ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)
ax2.spines["left"].set_linewidth(0.8)
ax2.spines["bottom"].set_linewidth(0.8)
ax2.tick_params(length=3, width=0.8, labelsize=9)

# Reference gene mean line for MCF10A
ref_line = dct_summary.loc["MCF10A", "mean"]
ax2.axhline(ref_line, color=COLORS["MCF10A"], linewidth=0.9,
            linestyle="--", zorder=2, alpha=0.7)
ax2.text(len(ORDER) - 0.45, ref_line + 0.05,
         "MCF10A\nreference", fontsize=7.5,
         color=COLORS["MCF10A"], style="italic")



plt.tight_layout()
plt.savefig("FigS2_delta_ct.png", dpi=300, bbox_inches="tight")
plt.savefig("FigS2_delta_ct.pdf", bbox_inches="tight")
print("Saved FigS2_delta_ct.png / .pdf")
plt.close()

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE S3 — Intra-assay replicate consistency heatmap
# ═══════════════════════════════════════════════════════════════════════════════

# Build matrix: rows = cell line + replicate, cols = target
hm_rows = []
for cell in ORDER:
    for rep in [1, 2, 3]:
        row = pivot[(pivot["Sample Name"] == cell) & (pivot["rep"] == rep)]
        if len(row) == 0:
            continue
        row = row.iloc[0]
        hm_rows.append({
            "Label":    f"{cell}\nrep {rep}",
            "Cell":     cell,
            "Rep":      rep,
            "SLC16A3":  row["SLC16A3"],
            "RNA18SN5": row["RNA18SN5"],
            "ACTB":     row["ACTB"],
        })

hm_df  = pd.DataFrame(hm_rows)

# Export raw Ct values to CSV
hm_df[["Cell", "Rep", "SLC16A3", "RNA18SN5", "ACTB"]].to_csv(
    "FigS3_raw_Ct_replicates.csv", index=False
)
print("Saved FigS3_raw_Ct_replicates.csv")

fig3, ax_sd = plt.subplots(1, 1, figsize=(6, 5.5))
fig3.patch.set_facecolor("white")

# ── SD per sample per target (reproducibility bars) ───────────────────────────
ax_sd.set_facecolor("white")

sd_rows = []
for cell in ORDER:
    sub = pivot[pivot["Sample Name"] == cell]
    for t in TARGETS:
        sd_rows.append({
            "Cell":   cell,
            "Target": t,
            "SD":     sub[t].std(ddof=1),
        })
sd_df = pd.DataFrame(sd_rows)

bar_h = 0.2
y_positions = []
ytick_labels = []

for ci, cell in enumerate(ORDER):
    base_y = ci * (len(TARGETS) * bar_h + 0.3)
    for ti, target in enumerate(TARGETS):
        y = base_y + ti * bar_h
        val = sd_df[(sd_df["Cell"] == cell) & (sd_df["Target"] == target)]["SD"].values[0]
        ax_sd.barh(y, val, height=bar_h * 0.75,
                   color=TARGET_COLORS[target],
                   edgecolor=EDGE, linewidth=0.6, alpha=0.88)
        ax_sd.text(val + 0.002, y, f"{val:.3f}",
                   va="center", fontsize=7.5, color=EDGE)
        y_positions.append(y)
        ytick_labels.append(f"$\\it{{{target}}}$")

# MIQE threshold line (SD ≤ 0.5 is acceptable)
ax_sd.axvline(0.5, color="#E74C3C", linewidth=1,
              linestyle="--", zorder=5)
ax_sd.text(0.51, ax_sd.get_ylim()[1] if ax_sd.get_ylim()[1] > 0 else 3.5,
           "MIQE\nthreshold\n(SD = 0.5)",
           fontsize=7, color="#E74C3C", va="top", style="italic")

ax_sd.set_yticks(y_positions)
ax_sd.set_yticklabels(ytick_labels, fontsize=8)
ax_sd.set_xlabel("Intra-assay SD of C$_t$", fontsize=9.5)
ax_sd.set_title("Replicate SD per cell line & target",
                fontsize=10.5, fontweight="bold", pad=8)
ax_sd.spines["top"].set_visible(False)
ax_sd.spines["right"].set_visible(False)
ax_sd.spines["left"].set_linewidth(0.8)
ax_sd.spines["bottom"].set_linewidth(0.8)
ax_sd.tick_params(length=3, width=0.8, labelsize=8.5)
ax_sd.set_xlim(0, 0.65)

# Cell line labels — placed above each group to avoid overlapping gene-name ticks
for ci, cell in enumerate(ORDER):
    base_y = ci * (len(TARGETS) * bar_h + 0.3)
    label_y = base_y + (len(TARGETS) - 1) * bar_h + 0.15
    ax_sd.text(-0.06, label_y, cell,
               transform=ax_sd.get_yaxis_transform(),
               ha="right", va="bottom",
               fontsize=8.5, fontstyle="italic",
               color=COLORS[cell], fontweight="bold")




plt.tight_layout()
plt.savefig("FigS3_replicate_heatmap.png", dpi=300, bbox_inches="tight")
plt.savefig("FigS3_replicate_heatmap.pdf", bbox_inches="tight")
print("Saved FigS3_replicate_heatmap.png / .pdf")
plt.close()


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE S4 — Melt temperature (Tm) distribution
# ═══════════════════════════════════════════════════════════════════════════════

import matplotlib.ticker as ticker

# Reload Tm1 data
df_tm = df[["Sample Name", "Target Name", "Cт"]].copy()  # df already loaded above

# Re-read with Tm1
df_raw2 = pd.read_excel(FILE, sheet_name="Results", engine="xlrd", header=None)
df2 = df_raw2.iloc[7:].copy()
df2.columns = df_raw2.iloc[7].values
df2 = df2.iloc[1:].reset_index(drop=True)
df2 = df2[df2["Well"].astype(str).str.match(r"^[A-Z]\d+$")].copy()

target_map2 = {"SLC": "SLC16A3", "18S": "RNA18SN5", "BACT": "ACTB"}
sample_map2 = {
    "MCF10A KONTROLA": "MCF10A",
    "MCF7":            "MCF7",
    "MDA231":          "MDA-MB-231",
    "MDA468":          "MDA-MB-468",
}
df2["Target Name"] = df2["Target Name"].map(target_map2)
df2["Sample Name"] = df2["Sample Name"].map(sample_map2)
df2["Tm1"] = pd.to_numeric(df2["Tm1"], errors="coerce")
df2 = df2[["Sample Name", "Target Name", "Tm1"]].dropna()

TARGETS   = ["SLC16A3", "RNA18SN5", "ACTB"]
ORDER     = ["MCF10A", "MCF7", "MDA-MB-231", "MDA-MB-468"]
COLORS    = {
    "MCF10A":    "#4A7FB5",
    "MCF7":      "#C0392B",
    "MDA-MB-231":"#E67E22",
    "MDA-MB-468":"#27AE60",
}
TARGET_COLORS = {
    "SLC16A3":  "#2C3E50",
    "RNA18SN5": "#8E44AD",
    "ACTB":     "#16A085",
}
EDGE = "#1a1a1a"
GREY = "#7f7f7f"
rng  = np.random.default_rng(42)

fig4, ax_dot = plt.subplots(1, 1, figsize=(7.5, 5))
fig4.patch.set_facecolor("white")
ax_dot.set_facecolor("white")

# Y positions: one group per target, spread cell lines within
target_y    = {"SLC16A3": 2, "RNA18SN5": 1, "ACTB": 0}
group_gap   = 0.18   # spread within target group

for target in TARGETS:
    base_y = target_y[target]
    overall_mean = df2[df2["Target Name"] == target]["Tm1"].mean()
    overall_std  = df2[df2["Target Name"] == target]["Tm1"].std(ddof=1)

    # Shaded band: mean ± SD
    ax_dot.axhspan(base_y - 0.38, base_y + 0.38,
                   color=TARGET_COLORS[target], alpha=0.07, zorder=1)

    # Mean line
    ax_dot.plot([overall_mean, overall_mean],
                [base_y - 0.32, base_y + 0.32],
                color=TARGET_COLORS[target], linewidth=1.8,
                zorder=3, alpha=0.6)

    for ci, cell in enumerate(ORDER):
        vals = df2[(df2["Target Name"] == target) &
                   (df2["Sample Name"] == cell)]["Tm1"].values
        y_offset = (ci - 1.5) * group_gap
        jitter   = rng.uniform(-0.03, 0.03, size=len(vals))

        ax_dot.scatter(vals, base_y + y_offset + jitter,
                       color=COLORS[cell], edgecolors=EDGE,
                       linewidths=0.7, s=55, zorder=4, alpha=0.92,
                       label=cell if target == "SLC16A3" else "")

    # Annotate mean ± SD
    ax_dot.text(overall_mean, base_y + 0.44,
                f"{overall_mean:.2f} ± {overall_std:.2f} °C",
                ha="center", va="bottom", fontsize=8.5,
                color=TARGET_COLORS[target], fontweight="bold")

# Y axis: target names
ax_dot.set_yticks([0, 1, 2])
ax_dot.set_yticklabels(
    [f"$\\it{{ACTB}}$\n(~{df2[df2['Target Name']=='ACTB']['Tm1'].mean():.1f} °C)",
     f"$\\it{{RNA18SN5}}$\n(~{df2[df2['Target Name']=='RNA18SN5']['Tm1'].mean():.1f} °C)",
     f"$\\it{{SLC16A3}}$\n(~{df2[df2['Target Name']=='SLC16A3']['Tm1'].mean():.1f} °C)"],
    fontsize=10
)
ax_dot.set_xlabel("Melting temperature T$_m$ (°C)", fontsize=10.5)
ax_dot.set_title("T$_m$ distribution per target and cell line",
                 fontsize=11, fontweight="bold", pad=10)
ax_dot.set_ylim(-0.6, 2.75)

# X axis: zoom to range of data with padding
x_min = df2["Tm1"].min() - 0.5
x_max = df2["Tm1"].max() + 0.5
ax_dot.set_xlim(x_min, x_max)
ax_dot.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))

ax_dot.spines["top"].set_visible(False)
ax_dot.spines["right"].set_visible(False)
ax_dot.spines["left"].set_linewidth(0.8)
ax_dot.spines["bottom"].set_linewidth(0.8)
ax_dot.tick_params(length=3, width=0.8)

# Separation between targets
for y in [0.5, 1.5]:
    ax_dot.axhline(y, color="#dddddd", linewidth=0.8, zorder=0)

# Legend
handles = [
    mpatches.Patch(facecolor=COLORS[c], edgecolor=EDGE, label=c)
    for c in ORDER
]
ax_dot.legend(handles=handles, fontsize=8.5, frameon=False,
              loc="upper center", bbox_to_anchor=(0.5, -0.18),
              ncol=4, handlelength=1.1,
              title="Cell line", title_fontsize=8.5)

# Export Tm summary statistics to CSV
rows = []
for target in TARGETS:
    for cell in ORDER:
        vals = df2[(df2["Target Name"] == target) &
                   (df2["Sample Name"] == cell)]["Tm1"].values
        rows.append({
            "Target":   target,
            "Cell line": cell,
            "Mean_Tm":  round(vals.mean(), 3),
            "SD":       round(vals.std(ddof=1), 3) if len(vals) > 1 else None,
            "Range":    f"{vals.min():.2f}–{vals.max():.2f}",
            "Peaks":    "1 (single)",
        })

pd.DataFrame(rows).to_csv("FigS4_Tm_summary.csv", index=False)
print("Saved FigS4_Tm_summary.csv")


plt.tight_layout(w_pad=2)
plt.savefig("FigS4_melt_temperature.png", dpi=300, bbox_inches="tight")
plt.savefig("FigS4_melt_temperature.pdf", bbox_inches="tight")
print("Saved FigS4_melt_temperature.png / .pdf")
plt.close()

print("\nAll four supplementary figures saved successfully.")
