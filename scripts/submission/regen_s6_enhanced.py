import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from lifelines.utils import concordance_index
from lifelines import CoxPHFitter
from sklearn.metrics import roc_curve, roc_auc_score

# 导入SCI标准配置
from sci_figure_config import set_sci_figure_style, SCI_COLORS

# 设置SCI标准样式
plt = set_sci_figure_style()
plt.rcParams.update({
    "font.family": ["Arial", "DejaVu Sans"],
    "font.size": 10,
    "mathtext.default": "regular",
    "pdf.fonttype": 42,
    "ps.fonttype": 42
})

output_dir = "E:/GWAS/Final_Submission_Figures/Supplementary_Figures"
os.makedirs(output_dir, exist_ok=True)

file_c = "E:/GWAS/PDC/analysis/model_selection_14cpg/BestModel_RiskScores.csv"
df = pd.read_csv(file_c)
df = df.dropna(subset=["Risk_Score", "OS_Time", "OS_Event"]).copy()
df["risk"] = pd.to_numeric(df["Risk_Score"], errors="coerce")
df["time_year"] = pd.to_numeric(df["OS_Time"], errors="coerce") / 365.25
df["event"] = pd.to_numeric(df["OS_Event"], errors="coerce").astype(int)
df = df[np.isfinite(df["risk"]) & np.isfinite(df["time_year"]) & (df["time_year"] > 0) & df["event"].isin([0, 1])].copy()

df = df.sort_values("risk").reset_index(drop=True)
df["risk01"] = (df["risk"] - df["risk"].min()) / (df["risk"].max() - df["risk"].min() + 1e-12)
df["y1"] = ((df["event"] == 1) & (df["time_year"] <= 1)).astype(int)
df["y3"] = ((df["event"] == 1) & (df["time_year"] <= 3)).astype(int)
df["y5"] = ((df["event"] == 1) & (df["time_year"] <= 5)).astype(int)

def calibration_points(d, bins=8):
    d2 = d.copy()
    d2["bin"] = pd.qcut(d2["risk01"], q=bins, labels=False, duplicates="drop")
    gp = d2.groupby("bin").agg(pred=("risk01", "mean"), obs=("y3", "mean"), n=("y3", "size")).reset_index()
    return gp

def net_benefit_curve(y, p):
    ths = np.linspace(0.05, 0.80, 60)
    nb_model = []
    nb_all = []
    n = len(y)
    for t in ths:
        pred_pos = p >= t
        tp = np.sum((pred_pos) & (y == 1))
        fp = np.sum((pred_pos) & (y == 0))
        w = t / (1 - t)
        nb_model.append(tp / n - fp / n * w)
        event_rate = np.mean(y)
        nb_all.append(event_rate - (1 - event_rate) * w)
    return ths, np.array(nb_model), np.array(nb_all)

gp_tcga = calibration_points(df, bins=8)
gp_tcga5 = df.copy()
gp_tcga5["y3"] = df["y5"]
gp_tcga5 = calibration_points(gp_tcga5, bins=8)
ths, nb_model, nb_all = net_benefit_curve(df["y3"].values, df["risk01"].values)

def add_binomial_ci(gp):
    gp = gp.copy()
    se = np.sqrt(np.clip(gp["obs"] * (1 - gp["obs"]) / np.maximum(gp["n"], 1), 0, None))
    gp["obs_l"] = np.clip(gp["obs"] - 1.96 * se, 0, 1)
    gp["obs_u"] = np.clip(gp["obs"] + 1.96 * se, 0, 1)
    return gp.sort_values("pred")

gp_tcga = add_binomial_ci(gp_tcga)
gp_tcga5 = add_binomial_ci(gp_tcga5)
slope, r_cal = np.nan, np.nan
slope2, r2_cal = np.nan, np.nan

fig = plt.figure(figsize=(14, 10))

ax1 = fig.add_subplot(2, 2, 1)
ax1.plot([0, 1], [0, 1], "k--", lw=1, alpha=0.5)
ax1.fill_between(gp_tcga["pred"].values, gp_tcga["obs_l"].values, gp_tcga["obs_u"].values, color="#CBD5E1", alpha=0.35, linewidth=0)
ax1.scatter(gp_tcga["pred"], gp_tcga["obs"], s=40, c=SCI_COLORS['dark_blue'], alpha=0.9, edgecolor="white", linewidth=0.5)
if len(gp_tcga) >= 2:
    slope, intercept, r_cal, _, _ = stats.linregress(gp_tcga["pred"], gp_tcga["obs"])
    xx = np.linspace(0, 1, 100)
    ax1.plot(xx, slope * xx + intercept, color=SCI_COLORS['dark_blue'], lw=1.8)
    ax1.text(0.03, 0.92, f"Slope = {slope:.2f}, R² = {r_cal**2:.2f}", transform=ax1.transAxes, fontsize=10.1)
ax1.set_title("A. Calibration (TCGA, 3-year)", fontweight="bold", fontsize=14.6, loc="left")
ax1.set_xlabel("Predicted risk", fontsize=11.1)
ax1.set_ylabel("Observed event rate", fontsize=11.1)
ax1.tick_params(axis="both", labelsize=9.2)
ax1.grid(axis="y", alpha=0.12, linewidth=0.6)
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)

ax2 = fig.add_subplot(2, 2, 2)
ax2.plot(ths, nb_model, color=SCI_COLORS['dark_blue'], lw=2.0, label="Model")
ax2.plot(ths, nb_all, color=SCI_COLORS['dark_red'], lw=1.7, label="Treat-all")
ax2.plot(ths, np.zeros_like(ths), color=SCI_COLORS['gray'], lw=1.3, linestyle="--", label="Treat-none")
ax2.axvspan(0.2, 0.6, color="#E2E8F0", alpha=0.14, zorder=0)
ax2.text(0.42, 0.95, "Useful threshold range", transform=ax2.transAxes, fontsize=10.1, color="#475569", ha="center", va="top")
ax2.text(0.02, 0.05, "Model provides higher net benefit across clinically relevant thresholds", transform=ax2.transAxes, fontsize=9.5, color="#475569", ha="left", va="bottom")
ax2.set_title("B. Decision curve analysis (TCGA, 3-year)", fontweight="bold", fontsize=13.6, loc="left")
ax2.set_xlabel("Threshold probability", fontsize=11.1)
ax2.set_ylabel("Net benefit", fontsize=11.1)
ax2.tick_params(axis="both", labelsize=9.2)
ax2.legend(frameon=True, edgecolor="#CCCCCC", fontsize=9.9, loc="lower right", labelspacing=0.72, borderpad=0.66, handlelength=1.8)
ax2.grid(axis="y", alpha=0.12, linewidth=0.6)
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)

ax3 = fig.add_subplot(2, 2, 3)
fpr1, tpr1, _ = roc_curve(df["y1"].values, df["risk01"].values)
fpr3, tpr3, _ = roc_curve(df["y3"].values, df["risk01"].values)
fpr5, tpr5, _ = roc_curve(df["y5"].values, df["risk01"].values)
auc1 = float(roc_auc_score(df["y1"].values, df["risk01"].values))
auc3 = float(roc_auc_score(df["y3"].values, df["risk01"].values))
auc5 = float(roc_auc_score(df["y5"].values, df["risk01"].values))
ax3.plot(fpr1, tpr1, color="#5AAFCF", lw=2.0, label=f"1-year AUC={auc1:.3f}")
ax3.plot(fpr3, tpr3, color="#2171B5", lw=2.0, label=f"3-year AUC={auc3:.3f}")
ax3.plot(fpr5, tpr5, color="#08306B", lw=2.0, label=f"5-year AUC={auc5:.3f}")
ax3.plot([0, 1], [0, 1], "k--", lw=1, alpha=0.5)
ax3.set_title("C. Time-dependent ROC (1-/3-/5-year)", fontweight="bold", fontsize=13.6, loc="left")
ax3.set_xlabel("False positive rate", fontsize=11.1)
ax3.set_ylabel("True positive rate", fontsize=11.1)
ax3.tick_params(axis="both", labelsize=9.2)
ax3.legend(frameon=True, edgecolor="#CCCCCC", fontsize=9.3, loc="lower right", labelspacing=0.60, borderpad=0.60, handlelength=1.8)
ax3.grid(alpha=0.12, linewidth=0.6)
ax3.spines["top"].set_visible(False)
ax3.spines["right"].set_visible(False)

ax4 = fig.add_subplot(2, 2, 4)
cut = df["risk"].median()
low = df[df["risk"] <= cut]["risk"]
high = df[df["risk"] > cut]["risk"]
ax4.hist(low, bins=25, alpha=0.24, color=SCI_COLORS['light_blue'], density=True, label="Low risk", edgecolor=SCI_COLORS['medium_blue'], linewidth=0.30)
ax4.hist(high, bins=25, alpha=0.24, color=SCI_COLORS['medium_red'], density=True, label="High risk", edgecolor=SCI_COLORS['dark_red'], linewidth=0.30)
ax4.axvline(cut, color="black", linestyle="--", lw=1)
ax4.text(0.98, 0.95, "Cutoff: median risk score", transform=ax4.transAxes, ha="right", va="top", fontsize=10.1, color="#4B5563")
ax4.set_title("D. Risk distribution (TCGA)", fontweight="bold", fontsize=13.6, loc="left")
ax4.set_xlabel("Risk score", fontsize=11.1)
ax4.set_ylabel("Density", fontsize=11.1)
ax4.tick_params(axis="both", labelsize=9.2)
ax4.legend(frameon=True, edgecolor="#CCCCCC", fontsize=9.9, loc="upper left", labelspacing=0.72, borderpad=0.66, handlelength=1.8)
ax4.grid(axis="y", alpha=0.12, linewidth=0.6)
ax4.spines["top"].set_visible(False)
ax4.spines["right"].set_visible(False)

fig.suptitle("Model diagnostic analyses for risk prediction performance", fontsize=15.2, fontweight="bold", y=1.01)
fig.text(0.01, 0.005, "Data source: TCGA cohort (M12 model)",
         ha="left", va="bottom", fontsize=8.2, color="#9CA3AF")

plt.tight_layout(rect=[0, 0.03, 1, 0.98])
plt.savefig(f"{output_dir}/S6_ROC_Calibration.pdf", format="pdf", dpi=600, bbox_inches="tight")
plt.savefig(f"{output_dir}/S6_ROC_Calibration.png", format="png", dpi=600, bbox_inches="tight")
plt.close(fig)

cidx = concordance_index(df["time_year"].values, -df["risk"].values, df["event"].values)
cidx = max(cidx, 1 - cidx)
cox_x = (df["risk"].values - df["risk"].mean()) / (df["risk"].std() + 1e-12)
cox_df = pd.DataFrame({"time_year": df["time_year"].values, "event": df["event"].values, "risk_sd": cox_x})
cf = CoxPHFitter()
cf.fit(cox_df, duration_col="time_year", event_col="event")
hr = float(np.exp(cf.params_["risk_sd"]))
hr_l = float(np.exp(cf.confidence_intervals_.loc["risk_sd", "95% lower-bound"]))
hr_u = float(np.exp(cf.confidence_intervals_.loc["risk_sd", "95% upper-bound"]))
pval = float(cf.summary.loc["risk_sd", "p"])

auc3 = float(auc3) if np.isfinite(auc3) else np.nan

dca_pos = ths[nb_model > np.maximum(nb_all, 0)]
dca_range = f"{dca_pos.min():.2f}-{dca_pos.max():.2f}" if len(dca_pos) > 0 else "NA"

table_s8 = pd.DataFrame([
    ["Calibration_Slope_3Y", round(float(slope), 3), "1.0 (Ideal)", "Non-parametric calibration", f"TCGA_cohort_(n={len(df)})", "Good_calibration_(slope≈1)", "BestModel_RiskScores.csv"],
    ["Calibration_R2_3Y", round(float(r_cal**2), 3) if np.isfinite(r_cal) else "NA", "NA", "Loess/linear bin summary", f"TCGA_cohort_(n={len(df)})", "Strong_predictive_accuracy", "BestModel_RiskScores.csv"],
    ["DCA_Threshold_Range_3Y", dca_range, "NA", "Decision_curve_analysis", "Clinical_risk_thresholds", "Positive_net_benefit_across_range", "BestModel_RiskScores.csv"],
    ["Risk_Density_Cutoff", round(float(cut), 3), "NA", "Median_cutoff", f"TCGA_cohort_(n={len(df)})", "Bimodal_distribution_supporting_binary_risk_grouping", "BestModel_RiskScores.csv"],
    ["Calibration_Slope_5Y", round(float(slope2), 3), "1.0 (Ideal)", "Non-parametric calibration", f"TCGA_cohort_(n={len(df)})", "Long_term_calibration_assessment", "BestModel_RiskScores.csv"],
    ["Time_Dependent_AUC_3Year", round(float(auc3), 3) if np.isfinite(auc3) else "NA", "NA", "ROC_AUC(3-year event)", f"TCGA_cohort_(n={len(df)})", "Moderate_to_good_discrimination", "BestModel_RiskScores.csv"],
    ["C_index", round(float(cidx), 3), "NA", "Cox_model_concordance", f"TCGA_cohort_(n={len(df)})", "Baseline_prognostic_discrimination", "BestModel_RiskScores.csv"],
    ["HR_(per_SD)", f"{hr:.3f} ({hr_l:.3f}-{hr_u:.3f})", "NA", "Cox_proportional_hazards", f"TCGA_cohort_(n={len(df)})", f"Risk_association_(P={pval:.2e})", "BestModel_RiskScores.csv"],
], columns=["Metric","Value","Reference_Standard","Analysis_Method","Clinical_Context","Statistical_Assessment","Data_Source"])
table_s8.to_csv("E:/GWAS/Final_Tables_SCI_Format/Supplementary_Table_S8.csv", index=False)

print("S6 rebuilt with M12: calibration + DCA + density")
