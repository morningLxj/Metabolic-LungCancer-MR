import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test
from scipy.stats import chi2
from sci_figure_config import set_sci_figure_style, SCI_COLORS

plt = set_sci_figure_style()
FS_PANEL = 9.0
FS_AXIS = 8.2
FS_TICK = 7.2
FS_TEXT = 8.0

coef_file = "E:/GWAS/PDC/analysis/model_selection_14cpg/BestModel_Coefficients.csv"
risk_file = "E:/GWAS/PDC/analysis/model_selection_14cpg/BestModel_RiskScores.csv"
cmp_file = "E:/GWAS/PDC/analysis/model_selection_14cpg/Model_Comparison_14CpG.csv"
out_dir = "E:/GWAS/Final_Submission_Figures/Main_Figures_14CpG"
os.makedirs(out_dir, exist_ok=True)

coef_df = pd.read_csv(coef_file)
risk_df = pd.read_csv(risk_file)
cmp_df = pd.read_csv(cmp_file)

risk_df["time_year"] = pd.to_numeric(risk_df["OS_Time"], errors="coerce") / 365.25
risk_df["event"] = pd.to_numeric(risk_df["OS_Event"], errors="coerce")
risk_df = risk_df[np.isfinite(risk_df["time_year"]) & np.isfinite(risk_df["event"]) & (risk_df["time_year"] > 0)]
risk_df["event"] = risk_df["event"].astype(int)
risk_df["Risk_Group"] = risk_df["Risk_Group"].astype(str)

cmp_df["N_CpG"] = pd.to_numeric(cmp_df["N_CpG"], errors="coerce")
cmp_df["AUC_Mean"] = pd.to_numeric(cmp_df["AUC_Mean"], errors="coerce")
cmp_df["C_Index_Oriented"] = pd.to_numeric(cmp_df["C_Index_Oriented"], errors="coerce")
cmp_df = cmp_df.sort_values("N_CpG")

cox_data = risk_df[["time_year", "event", "Risk_Group"]].copy()
cox_data["Risk_Group_Binary"] = (cox_data["Risk_Group"] == "High Risk").astype(int)
cph = CoxPHFitter()
cph.fit(cox_data[["time_year", "event", "Risk_Group_Binary"]], duration_col="time_year", event_col="event")
hr = float(np.exp(cph.params_["Risk_Group_Binary"]))
hr_l = float(np.exp(cph.confidence_intervals_.loc["Risk_Group_Binary", "95% lower-bound"]))
hr_u = float(np.exp(cph.confidence_intervals_.loc["Risk_Group_Binary", "95% upper-bound"]))

g_high = risk_df[risk_df["Risk_Group"] == "High Risk"]
g_low = risk_df[risk_df["Risk_Group"] == "Low Risk"]
lr = logrank_test(
    g_high["time_year"], g_low["time_year"],
    event_observed_A=g_high["event"], event_observed_B=g_low["event"]
)
km_p = float(lr.p_value)

fig, axes = plt.subplots(2, 2, figsize=(11.6, 9.2))
axA, axB, axC, axD = axes.flatten()

for ax in [axA, axB, axC, axD]:
    ax.set_facecolor("white")

axA.set_title("A. Causal-guided feature selection pipeline", loc="left", fontweight="bold", fontsize=FS_PANEL)
axA.axis("off")
nodes = [
    (0.50, 0.85, "GWAS / MR /\nColocalization", "#F8FAFC", "#1F2937"),
    (0.50, 0.64, "Candidate\nCpGs (n=14)", "#EFF6FF", "#1E3A8A"),
    (0.50, 0.43, "Penalized Cox\nselection", "#FFF7ED", "#9A3412"),
    (0.50, 0.22, "Final model\nM12 (12 CpGs)", "#EEF2FF", "#312E81"),
]
for x, y, txt, fc, tc in nodes:
    box = FancyBboxPatch((x - 0.15, y - 0.085), 0.30, 0.165, boxstyle="round,pad=0.02", fc=fc, ec="#94A3B8", lw=1.0, transform=axA.transAxes)
    axA.add_patch(box)
    axA.text(x, y, txt, ha="center", va="center", fontsize=FS_TEXT, color=tc, transform=axA.transAxes)
for y0, y1 in [(0.76, 0.71), (0.55, 0.50), (0.34, 0.29)]:
    axA.annotate("", xy=(0.50, y1), xytext=(0.50, y0), xycoords=axA.transAxes, textcoords=axA.transAxes, arrowprops=dict(arrowstyle="-|>", color="#64748B", lw=1.0))

axB.set_title("B. Model-size selection profile", loc="left", fontweight="bold", fontsize=FS_PANEL)
axB.plot(cmp_df["N_CpG"], cmp_df["AUC_Mean"], color="#2166AC", marker="o", lw=1.2, label="Mean AUC")
axB.plot(cmp_df["N_CpG"], cmp_df["C_Index_Oriented"], color="#B2182B", marker="^", lw=1.2, linestyle="--", label="C-index")
axB.axvline(12, color="#B2182B", linestyle=":", lw=1.2)
axB.text(
    12.12,
    float(np.nanmax(np.r_[cmp_df["AUC_Mean"].values, cmp_df["C_Index_Oriented"].values])) + 0.002,
    "M12 selected",
    color="#B2182B",
    fontsize=FS_TICK,
    bbox=dict(boxstyle="round,pad=0.18", facecolor="white", edgecolor="none", alpha=0.88),
)
axB.set_xlabel("Number of CpGs", fontsize=FS_AXIS)
axB.set_ylabel("Performance", fontsize=FS_AXIS)
axB.tick_params(axis="both", labelsize=FS_TICK)
axB.grid(alpha=0.18, axis="y")
axB.legend(frameon=False, fontsize=FS_TICK, loc="lower right")

coef_df["abs_coef"] = coef_df["Coef"].abs()
coef_df = coef_df.sort_values("abs_coef", ascending=True)
bar_colors = np.where(coef_df["Coef"] >= 0, SCI_COLORS["dark_red"], SCI_COLORS["medium_blue"])
axC.barh(coef_df["Feature"], coef_df["Coef"], color=bar_colors, alpha=0.92)
axC.axvline(0, color="gray", linestyle="--", lw=1.0)
axC.set_title("C. Regression coefficients of the final 12-CpG model", loc="left", fontweight="bold", fontsize=FS_PANEL)
axC.set_xlabel("Coefficient (β)", fontsize=FS_AXIS)
axC.set_ylabel("CpG", fontsize=FS_AXIS)
axC.grid(alpha=0.18, axis="x")
axC.tick_params(axis="x", labelsize=FS_TICK)
axC.tick_params(axis="y", labelsize=6.8)
axC.margins(y=0.05)

kmh = KaplanMeierFitter()
kml = KaplanMeierFitter()
kmh.fit(g_high["time_year"], event_observed=g_high["event"], label=f"High Risk (n={len(g_high)})")
kml.fit(g_low["time_year"], event_observed=g_low["event"], label=f"Low Risk (n={len(g_low)})")
kmh.plot_survival_function(ax=axD, ci_show=False, color="#B2182B", lw=1.4)
kml.plot_survival_function(ax=axD, ci_show=False, color="#2166AC", lw=1.4)
axD.set_title("D. Risk stratification in the training cohort", loc="left", fontweight="bold", fontsize=FS_PANEL)
axD.set_xlabel("Time (Years)", fontsize=FS_AXIS)
axD.set_ylabel("Overall Survival Probability", fontsize=FS_AXIS)
axD.tick_params(axis="both", labelsize=FS_TICK)
axD.grid(alpha=0.18, axis="y")
axD.legend(frameon=False, fontsize=FS_TICK, loc="upper right")
axD.text(
    risk_df["time_year"].max() * 0.955, 0.84,
    f"HR={hr:.2f} ({hr_l:.2f}-{hr_u:.2f})",
    ha="right", va="top", fontsize=FS_TEXT,
    bbox=dict(boxstyle="round,pad=0.22", facecolor="white", edgecolor="#CBD5E1")
)
axD.text(
    risk_df["time_year"].max() * 0.955, 0.75,
    f"Log-rank P={km_p:.2e}",
    ha="right", va="top", fontsize=FS_TEXT,
    bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="#CBD5E1")
)

fig.suptitle("Model construction and training-cohort stratification (M12)", fontsize=11.6, fontweight="bold", y=0.988)
fig.subplots_adjust(left=0.07, right=0.99, top=0.93, bottom=0.07, wspace=0.25, hspace=0.32)

fig.savefig(f"{out_dir}/Figure5_M12_Standalone_py.tif", format="tiff", dpi=600)
fig.savefig(f"{out_dir}/Figure5_M12_Standalone_py.png", format="png", dpi=600)
fig.savefig(f"{out_dir}/Figure5_M12_Standalone_py.pdf", format="pdf", dpi=600)
plt.close(fig)

pd.DataFrame({
    "Model": ["M12_Top12"],
    "N": [len(risk_df)],
    "Events": [int(risk_df["event"].sum())],
    "HR_High_vs_Low": [hr],
    "HR_CI_Lower": [hr_l],
    "HR_CI_Upper": [hr_u],
    "KM_P": [km_p]
}).to_csv(f"{out_dir}/Figure5_M12_Summary_py.csv", index=False)

print("Figure5 M12 python standalone generated")
