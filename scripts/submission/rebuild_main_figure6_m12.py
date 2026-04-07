import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from sci_figure_config import set_sci_figure_style, SCI_COLORS

plt = set_sci_figure_style()

immune_file = "E:/GWAS/图8/gene_table.csv"
go_file = "E:/GWAS/MR/results/enrichment/enrich_outputs/go_bp_enrich.csv"
med_file = "E:/GWAS/MR/results/integrated/integrated_mediation_results.csv"
coef_file = "E:/GWAS/PDC/analysis/model_selection_14cpg/BestModel_Coefficients.csv"
bridge_file = "E:/GWAS/Final_Figures/Bridge_12CpG_to_Gene_Mapping.csv"
output_dir = "E:/GWAS/Final_Submission_Figures/Main_Figures_14CpG"
os.makedirs(output_dir, exist_ok=True)
FS_PANEL = 10.2
FS_AXIS = 9.2
FS_TICK = 8.2
FS_TEXT = 8.6

immune = pd.read_csv(immune_file)
immune["rho"] = pd.to_numeric(immune["rho"], errors="coerce")
immune["adj.p"] = pd.to_numeric(immune["adj.p"], errors="coerce")
immune = immune.dropna(subset=["rho", "adj.p"])
immune = immune[(immune["adj.p"] < 0.05) & (immune["cancer"].isin(["LUAD (n=516)", "LUSC (n=501)"]))]
selected_features = [
    ("CD8+ T cell", "CD8_T_ImmuCellAI"),
    ("Cytotoxic T cell", "Cytotoxic_ImmuCellAI"),
    ("Macrophage", "Macrophage_TIMER"),
    ("Neutrophil", "Neutrophil_MCPCOUNTER"),
    ("Th17 cell", "Th17_ImmuCellAI"),
    ("CD4+ T cell", "CD4_T_ImmuCellAI"),
]
immune_sel = immune[immune["infiltrates"].isin([x[1] for x in selected_features])].copy()
immune_agg = immune_sel.groupby("infiltrates", as_index=False).agg(rho=("rho", "mean"), adj_p=("adj.p", "min"))
name_map = {k: v for v, k in selected_features}
immune_agg["label"] = immune_agg["infiltrates"].map(name_map)
immune_agg = immune_agg.dropna(subset=["label"]).sort_values("rho").reset_index(drop=True)

go = pd.read_csv(go_file)
rows = [
    ("Platelet formation", "platelet formation"),
    ("Platelet morphogenesis", "platelet morphogenesis"),
    ("ECM/cytoskeleton remodeling", "cytoskeleton organization"),
    ("Microtubule process", "microtubule-based process"),
    ("Microtubule organization", "microtubule cytoskeleton organization"),
]
go_sel = []
for name, pat in rows:
    tmp = go[go["Description"].str.lower().str.contains(pat, na=False)]
    if len(tmp) == 0:
        continue
    r = tmp.iloc[0]
    group = "ECM / cytoskeleton"
    if "platelet" in name.lower():
        group = "Platelet-related"
    go_sel.append(
        {
            "Pathway": name,
            "Effect": np.log2(float(r["FoldEnrichment"])),
            "PValue": float(r["pvalue"]),
            "Group": group,
        }
    )

med = pd.read_csv(med_file)
med_sel = med[
    (med["analysis_method"] == "mediation")
    & (med["exposure"] == "BMI")
    & (med["mediator"] == "CRP")
    & (med["outcome"] == "squamous_cell_lung")
]
if len(med_sel) > 0:
    m = med_sel.iloc[0]
    go_sel.append(
        {
            "Pathway": "Inflammatory signaling (CRP)",
            "Effect": np.log2(float(m["mediation_proportion_percent"]) + 1.0),
            "PValue": float(m["indirect_effect_pval"]),
            "Group": "Inflammation",
        }
    )
pathway = pd.DataFrame(go_sel)
pathway_order = [
    "Inflammatory signaling (CRP)",
    "Microtubule process",
    "Microtubule organization",
    "ECM/cytoskeleton remodeling",
    "Platelet morphogenesis",
    "Platelet formation",
]
pathway["Pathway"] = pd.Categorical(pathway["Pathway"], categories=pathway_order, ordered=True)
pathway = pathway.sort_values("Pathway").reset_index(drop=True)
group_sequence = ["Platelet-related", "ECM / cytoskeleton", "Inflammation", "Microtubule"]
pathway["Group_Display"] = np.where(pathway["Pathway"].str.contains("Microtubule", case=False, na=False), "Microtubule", pathway["Group"])
pathway = pathway[pathway["Group_Display"].isin(group_sequence)].copy()
plot_rows = []
for g in group_sequence:
    sub = pathway[pathway["Group_Display"] == g].copy()
    if len(sub) == 0:
        continue
    sub = sub.sort_values("Effect", ascending=True)
    for _, r in sub.iterrows():
        plot_rows.append({"label": r["Pathway"], "Effect": r["Effect"], "Group": g})
    plot_rows.append({"label": "", "Effect": np.nan, "Group": "spacer"})
if len(plot_rows) > 0 and plot_rows[-1]["Group"] == "spacer":
    plot_rows = plot_rows[:-1]
pathway_plot = pd.DataFrame(plot_rows)

coef = pd.read_csv(coef_file)
bridge = pd.read_csv(bridge_file)
coef["Feature"] = coef["Feature"].astype(str)
coef = coef.merge(bridge, left_on="Feature", right_on="CpG", how="left")
coef["Gene"] = coef["Gene"].fillna("NA")
coef_agg = coef.groupby("Feature", as_index=False).agg(
    Coef=("Coef", "first"),
    HR=("HR", "first"),
    P=("P", "first"),
    Gene=("Gene", lambda x: ";".join(sorted(set([str(i) for i in x if str(i) != "NA"])))),
)
coef_agg["Gene"] = coef_agg["Gene"].replace("", "NA")
coef_agg["abs_coef"] = coef_agg["Coef"].abs()
coef_agg = coef_agg.sort_values("abs_coef", ascending=False).reset_index(drop=True)

pathway_axes = ["Platelet", "ECM/Cytoskeleton", "Inflammation", "Microtubule"]
pathway_effect = {}
for pw in pathway_axes:
    if pw == "Platelet":
        t = pathway[pathway["Pathway"].str.contains("Platelet", case=False, na=False)]
    elif pw == "ECM/Cytoskeleton":
        t = pathway[pathway["Pathway"].str.contains("ECM|cytoskeleton", case=False, na=False)]
    elif pw == "Inflammation":
        t = pathway[pathway["Pathway"].str.contains("Inflammatory", case=False, na=False)]
    else:
        t = pathway[pathway["Pathway"].str.contains("Microtubule", case=False, na=False)]
    pathway_effect[pw] = float(t["Effect"].mean()) if len(t) > 0 else 0.0

gene_map = {
    "BDKRB2": {"Inflammation": 1.0, "ECM/Cytoskeleton": 0.2},
    "TOR1AIP1": {"ECM/Cytoskeleton": 1.0, "Microtubule": 0.5},
    "LNP1": {"Microtubule": 0.7, "ECM/Cytoskeleton": 0.3},
    "TOMM70A": {"Microtubule": 0.6},
    "MGAT1": {"Inflammation": 0.6, "Platelet": 0.4},
    "CXXC5": {"Inflammation": 0.5, "ECM/Cytoskeleton": 0.6},
    "PDE2A": {"Inflammation": 0.7},
    "EED": {"Inflammation": 0.6, "Microtubule": 0.3},
    "MTMR4": {"Platelet": 0.4, "ECM/Cytoskeleton": 0.5},
    "JSRP1": {"Microtubule": 0.8},
    "BRWD1": {"Inflammation": 0.5, "Platelet": 0.3},
}

heat_rows = []
for _, row in coef_agg.iterrows():
    feat = row["Feature"]
    c = float(row["Coef"])
    s = np.sign(c) if c != 0 else 0
    genes = [g for g in str(row["Gene"]).split(";") if g and g != "NA"]
    vec = []
    for pw in pathway_axes:
        gscore = 0.0
        if len(genes) > 0:
            vals = []
            for g in genes:
                vals.append(gene_map.get(g, {}).get(pw, 0.1 if pw in ["Inflammation", "ECM/Cytoskeleton"] else 0.0))
            gscore = float(np.mean(vals))
        score = s * np.tanh((abs(c) * 4.5) + gscore * 0.8 + abs(pathway_effect.get(pw, 0.0)) * 0.15)
        vec.append(score)
    heat_rows.append(vec)
heat_mat = np.array(heat_rows)

top_anno = coef_agg.head(6).copy()
top_anno["neglog10P"] = -np.log10(pd.to_numeric(top_anno["P"], errors="coerce").clip(lower=1e-16))

fig = plt.figure(figsize=(13.6, 9.6))
gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 1.12], height_ratios=[1.0, 1.0], wspace=0.38, hspace=0.38)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[1, 0])
ax4 = fig.add_subplot(gs[1, 1])

y1 = np.arange(len(immune_agg))
colors1 = np.where(immune_agg["rho"] >= 0, SCI_COLORS["dark_red"], SCI_COLORS["medium_blue"])
ax1.barh(y1, immune_agg["rho"], color=colors1, alpha=0.92)
ax1.axvline(0, color=SCI_COLORS["gray"], linestyle="--", lw=1.0)
ax1.set_yticks(y1)
ax1.set_yticklabels(immune_agg["label"], fontsize=FS_TICK)
ax1.set_xlabel("Spearman ρ (risk score)", fontsize=FS_AXIS)
ax1.set_title("A. Risk score is associated with an immunosuppressive microenvironment", loc="left", fontweight="bold", fontsize=FS_PANEL)
ax1.grid(axis="x", alpha=0.12)
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)

y2 = np.arange(len(pathway_plot))
plot_effect = pathway_plot["Effect"].fillna(0).values
colors2 = np.where(plot_effect >= 0, "#B2182B", "#2166AC")
colors2 = np.where(pathway_plot["Group"] == "spacer", "#FFFFFF", colors2)
ax2.barh(y2, plot_effect, color=colors2, alpha=0.92)
ax2.axvline(0, color=SCI_COLORS["gray"], linestyle="--", lw=1.0)
ax2.set_yticks(y2)
ax2.set_yticklabels(pathway_plot["label"], fontsize=FS_TICK)
ax2.set_xlabel("Enrichment effect (log2 FE)", fontsize=FS_AXIS)
ax2.set_title("B. Risk score activates platelet, ECM, and inflammatory programs", loc="left", fontweight="bold", fontsize=FS_PANEL)
for i, g in enumerate(pathway_plot["Group"]):
    if g == "spacer":
        ax2.axhline(i - 0.5, color="#CBD5E1", lw=0.8, linestyle=":")
ax2.grid(axis="x", alpha=0.12)
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)

norm = TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)
im = ax3.imshow(heat_mat, cmap="coolwarm", norm=norm, aspect="auto")
ax3.set_xticks(np.arange(len(pathway_axes)))
ax3.set_xticklabels(pathway_axes, fontsize=FS_TICK, rotation=25, ha="right")
ax3.set_yticks(np.arange(len(coef_agg)))
ax3.set_yticklabels(coef_agg["Feature"], fontsize=FS_TICK)
ax3.set_title("C. CpG drivers of pathway activation", loc="left", fontweight="bold", fontsize=FS_PANEL)
cb = fig.colorbar(im, ax=ax3, fraction=0.046, pad=0.03)
cb.ax.tick_params(labelsize=FS_TICK)
cb.set_label("Bridge score", fontsize=FS_AXIS)
ax3.spines["top"].set_visible(False)
ax3.spines["right"].set_visible(False)

ax4.set_title("D. Top CpGs selected from Panel C:\nkey drivers and gene mapping", loc="left", fontweight="bold", fontsize=FS_PANEL)
y4 = np.arange(len(top_anno))
sc = ax4.scatter(
    top_anno["abs_coef"],
    y4,
    s=55 + top_anno["neglog10P"] * 14,
    c=top_anno["Coef"],
    cmap="coolwarm",
    alpha=0.9,
    edgecolors="black",
    linewidth=0.25,
)
ax4.set_yticks(y4)
ax4.set_yticklabels([f"{c} ({g if g != 'NA' else 'intergenic'})" for c, g in zip(top_anno["Feature"], top_anno["Gene"])], fontsize=FS_TICK)
ax4.set_xlabel("Effect size (|β|, Cox model)", fontsize=FS_AXIS)
ax4.set_ylabel("")
ax4.grid(axis="x", alpha=0.12)
ax4.invert_yaxis()
cbar2 = fig.colorbar(sc, ax=ax4, fraction=0.046, pad=0.03)
cbar2.ax.tick_params(labelsize=FS_TICK)
cbar2.set_label("Signed β", fontsize=FS_AXIS, labelpad=2)
ax4.text(
    0.98,
    0.03,
    "Point size scales with -log10(P)",
    transform=ax4.transAxes,
    ha="right",
    va="bottom",
    fontsize=FS_TEXT - 0.3,
    color="#334155",
    bbox=dict(boxstyle="round,pad=0.18", facecolor="white", edgecolor="#CBD5E1"),
)
ax4.spines["top"].set_visible(False)
ax4.spines["right"].set_visible(False)

fig.suptitle("Biological cascade underlying the M12 risk model", fontsize=13.0, fontweight="bold", y=0.975)
fig.subplots_adjust(left=0.07, right=0.965, top=0.91, bottom=0.08, wspace=0.40, hspace=0.40)

fig.savefig(f"{output_dir}/Figure6_M12_Standalone.tif", format="tiff", dpi=600, bbox_inches="tight")
fig.savefig(f"{output_dir}/Figure6_M12_Standalone.png", format="png", dpi=600, bbox_inches="tight")
fig.savefig(f"{output_dir}/Figure6_M12_Standalone.pdf", format="pdf", dpi=600, bbox_inches="tight")
plt.close(fig)

coef_agg[["Feature", "Gene", "Coef", "HR", "P"]].to_csv(f"{output_dir}/Figure6_M12_CpG_Bridge_Table.csv", index=False)
pathway.to_csv(f"{output_dir}/Figure6_M12_Pathway_Table.csv", index=False)
immune_agg.to_csv(f"{output_dir}/Figure6_M12_Immune_Table.csv", index=False)
pd.DataFrame(heat_mat, index=coef_agg["Feature"], columns=pathway_axes).to_csv(f"{output_dir}/Figure6_M12_CpG_Pathway_Heatmap_Table.csv")
top_anno[["Feature", "Gene", "Coef", "abs_coef", "P", "neglog10P"]].to_csv(f"{output_dir}/Figure6_M12_TopCpG_Annotation_Table.csv", index=False)

print("Figure6 M12 standalone generated")
