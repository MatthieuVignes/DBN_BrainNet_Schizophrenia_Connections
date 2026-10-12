"""Figures for the results summary, from the outputs of scripts 24/25.
Usage: python3 27_make_summary_figures.py [DIR_WITH_CSVS] [OUT_DIR]"""
import sys, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
d = sys.argv[1] if len(sys.argv) > 1 else "."
o = sys.argv[2] if len(sys.argv) > 2 else "."
full = pd.read_csv(f"{d}/bootstrap_difference_test_results.csv")
ns = pd.read_csv(f"{d}/bootstrap_difference_nonsmoker.csv")
cnt = pd.read_csv(f"{d}/per_bootstrap_edge_counts.csv")
C = {"HC": "#4C78A8", "PDS": "#E45756"}
# Fig 2: distribution of edge counts across bootstrap fits
fig, ax = plt.subplots(figsize=(6.4, 3.8))
groups = ["H.E", "S.E", "H.R", "S.R"]; pos = [1, 2, 3.4, 4.4]
bp = ax.boxplot([cnt[g].dropna() for g in groups], positions=pos, widths=.6, patch_artist=True, showfliers=False)
for patch, g in zip(bp["boxes"], groups):
    patch.set_facecolor(C["HC"] if g[0] == "H" else C["PDS"]); patch.set_alpha(.75)
for m in bp["medians"]: m.set_color("black")
ax.set_xticks([1.5, 3.9]); ax.set_xticklabels(["Encoding", "Retrieval"])
ax.set_ylabel("Edges per bootstrap fit"); ax.set_title("Edge count across 100 bootstrap fits per group", fontsize=11)
ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, fc=C["HC"], alpha=.75), plt.Rectangle((0, 0), 1, 1, fc=C["PDS"], alpha=.75)],
          labels=["HC", "PDS"], frameon=False, loc="upper center", ncol=2, bbox_to_anchor=(0.5, -0.12))
for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig(f"{o}/fig2_edge_count_distribution.png", dpi=200); plt.close(fig)
# Fig 3: forest plot of contrasts (full sample and non-smokers)
fig, ax = plt.subplots(figsize=(7.2, 3.6))
labels = ["Encoding: PDS − HC", "Retrieval: PDS − HC", "HC: Retrieval − Encoding", "PDS: Retrieval − Encoding"]
for i, lab in enumerate(labels):
    for df, off, col, name in ((full, .15, "#333333", "All participants"), (ns, -.15, "#F58518", "Non-smokers")):
        r = df.iloc[i]; y = len(labels) - 1 - i + off
        ax.plot([r.ci_lower, r.ci_upper], [y, y], color=col, lw=2)
        ax.plot(r.mean_diff_edges, y, "o", color=col, label=name if i == 0 else None)
ax.axvline(0, color="grey", ls="--", lw=1)
ax.set_yticks(range(len(labels))); ax.set_yticklabels(labels[::-1])
ax.set_xlabel("Difference in edge count (95% bootstrap interval)")
ax.set_title("No contrast excludes zero", fontsize=11)
ax.legend(frameon=False, loc="upper center", ncol=2, bbox_to_anchor=(0.4, -0.2))
for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig(f"{o}/fig3_forest_plot.png", dpi=200); plt.close(fig)
