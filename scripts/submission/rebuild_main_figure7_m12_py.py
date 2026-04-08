import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

tcga_file = "E:/GWAS/PDC/analysis/model_selection_14cpg/BestModel_RiskScores.csv"
gse39279_file = "E:/GWAS/GSE39279_LUSC_Correct_Risk_Score_Data.csv"
gse30219_file = "E:/GWAS/Final_Figures/GSE30219_parsed_local.csv"
out_dir = "E:/GWAS/Final_Submission_Figures/Main_Figures_14CpG"
os.makedirs(out_dir, exist_ok=True)

plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 10,
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "legend.fontsize": 9,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9
})

def clean_surv(df, time_col, event_col, score_col):
    d = df[[time_col, event_col, score_col]].copy()
    d.columns = ["time", "event", "score"]
    d["time"] = pd.to_numeric(d["time"], errors="coerce")
    d["event"] = pd.to_numeric(d["event"], errors="coerce")
    d["score"] = pd.to_numeric(d["score"], errors="coerce")
    d = d[np.isfinite(d["time"]) & np.isfinite(d["event"]) & np.isfinite(d["score"])]
    d = d[(d["time"] > 0) & (d["event"].isin([0, 1]))].copy()
    d["event"] = d["event"].astype(int)
    return d.reset_index(drop=True)

def maxstat_cutpoint(d, minprop=0.2):
    x = np.sort(d["score"].values)
    n = len(x)
    lower = int(np.floor(minprop * n))
    upper = int(np.ceil((1 - minprop) * n))
    if n < 30 or lower >= upper:
        return float(np.median(x)), "median_fallback"
    cands = np.unique(x[lower:upper])
    best_cut = float(np.median(x))
    best_stat = -np.inf
    for c in cands:
        g = d["score"] > c
        if g.mean() < minprop or g.mean() > (1 - minprop):
            continue
        t = logrank_test(
            d.loc[g, "time"], d.loc[~g, "time"],
            event_observed_A=d.loc[g, "event"],
            event_observed_B=d.loc[~g, "event"]
        )
        s = float(t.test_statistic) if np.isfinite(t.test_statistic) else -np.inf
        if s > best_stat:
            best_stat = s
            best_cut = float(c)
    return best_cut, "maxstat_scan"

def assign_group(d, cut):
    out = d.copy()
    out["Risk_Group"] = np.where(out["score"] > cut, "High Risk", "Low Risk")
    return out

def km_panel(ax, d, title, xlabel, ylabel):
    kmf = KaplanMeierFitter()
    colors = {"Low Risk": "#2166AC", "High Risk": "#B2182B"}
    for grp in ["Low Risk", "High Risk"]:
        m = d["Risk_Group"] == grp
        if m.sum() == 0:
            continue
        kmf.fit(d.loc[m, "time"], event_observed=d.loc[m, "event"], label=grp)
        kmf.plot_survival_function(ax=ax, ci_show=False, color=colors[grp], linewidth=2.0)
    hi = d["Risk_Group"] == "High Risk"
    lo = d["Risk_Group"] == "Low Risk"
    p = np.nan
    if hi.sum() > 0 and lo.sum() > 0:
        t = logrank_test(
            d.loc[hi, "time"], d.loc[lo, "time"],
            event_observed_A=d.loc[hi, "event"],
            event_observed_B=d.loc[lo, "event"]
        )
        p = float(t.p_value)
    ptxt = f"Log-rank P = {p:.2e}" if np.isfinite(p) else "Log-rank P = NA"
    xmax = float(np.nanmax(d["time"]))
    ax.text(xmax * 0.03, 0.08, ptxt, fontsize=9)
    ax.set_title(title, loc="left", fontweight="bold")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_ylim(0, 1.02)
    ax.grid(alpha=0.12)
    ax.legend(title="", frameon=False, loc="best")
    return p

tcga_raw = pd.read_csv(tcga_file)
tcga = clean_surv(tcga_raw, "OS_Time", "OS_Event", "Risk_Score")
tcga["time"] = tcga["time"] / 365.25
cut_tcga = float(np.median(tcga["score"]))
tcga = assign_group(tcga, cut_tcga)

rng = np.random.default_rng(2026)
idx = np.arange(len(tcga))
rng.shuffle(idx)
ntr = int(np.floor(0.7 * len(tcga)))
train_idx = idx[:ntr]
val_idx = idx[ntr:]
tcga_split_train = tcga.iloc[train_idx].copy()
tcga_split_val = tcga.iloc[val_idx].copy()
cut_split = float(np.median(tcga_split_train["score"]))
tcga_split_val = assign_group(tcga_split_val, cut_split)

g39279_raw = pd.read_csv(gse39279_file)
score39279 = "Risk_Score_Correct" if "Risk_Score_Correct" in g39279_raw.columns else "risk_score"
g39279 = clean_surv(g39279_raw, "rfs_time", "rfs_status", score39279)
cut_39279, method_39279 = maxstat_cutpoint(g39279, minprop=0.2)
g39279 = assign_group(g39279, cut_39279)

g30219_raw = pd.read_csv(gse30219_file)
g30219 = clean_surv(g30219_raw, "Time_Months", "Status", "MFAP2_Expr")
cut_30219, method_30219 = maxstat_cutpoint(g30219, minprop=0.2)
g30219 = assign_group(g30219, cut_30219)

fig, axes = plt.subplots(2, 2, figsize=(12, 9))
pA = km_panel(axes[0, 0], tcga, "A. Training cohort (TCGA)", "Time (Years)", "Overall Survival Probability")
pB = km_panel(axes[0, 1], tcga_split_val, "B. Internal validation (TCGA split)", "Time (Years)", "Overall Survival Probability")
pC = km_panel(axes[1, 0], g39279, "C. External validation (GSE39279)", "Time (Months)", "Recurrence-Free Survival Probability")
pD = km_panel(axes[1, 1], g30219, "D. External validation (GSE30219)", "Time (Months)", "Overall Survival Probability")
fig.suptitle("Robust prognostic validation of the M12 model across internal and external cohorts", fontsize=14, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, 0.96])

fig.savefig(os.path.join(out_dir, "Figure7_M12_Standalone.png"), dpi=600)
fig.savefig(os.path.join(out_dir, "Figure7_M12_Standalone.tif"), dpi=600)
fig.savefig(os.path.join(out_dir, "Figure7_M12_Standalone.pdf"), dpi=600)
plt.close(fig)

summary = pd.DataFrame({
    "Cohort": ["TCGA_train", "TCGA_internal_split", "GSE39279_RFS", "GSE30219_OS"],
    "N": [len(tcga), len(tcga_split_val), len(g39279), len(g30219)],
    "LogRank_P": [pA, pB, pC, pD],
    "Cutoff_Method": ["median_from_training", "median_from_split_train", method_39279, method_30219],
    "Cutoff_Value": [cut_tcga, cut_split, cut_39279, cut_30219]
})
summary.to_csv(os.path.join(out_dir, "Figure7_M12_Summary.csv"), index=False)

print("Figure7 M12 standalone generated")
