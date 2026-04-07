import os
import re
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from sklearn.metrics import roc_auc_score

base = "E:/GWAS"
out_dir = f"{base}/Final_Tables_SCI_Format"
os.makedirs(out_dir, exist_ok=True)

# Inputs
s1_src_file = f"{base}/MK2025112604E/投稿/MRpaper1/Supplementary_Tables_From_PDF/Supplementary_Table_S1.csv"
s2_src_file = f"{base}/MK2025112604E/投稿/MRpaper1/Supplementary_Tables_From_PDF/Supplementary_Table_S2.csv"
s3_src_file = f"{base}/MK2025112604E/投稿/MRpaper1/Supplementary_Tables_From_PDF/Supplementary_Table_S3.csv"
coloc_file = f"{base}/Table_S4_Enhanced_Colocalization_Analysis.csv"
bridge_file = f"{base}/Final_Figures/Bridge_12CpG_to_Gene_Mapping.csv"
coef_file = f"{base}/PDC/analysis/model_selection_14cpg/BestModel_Coefficients.csv"
lasso_sum_file = f"{base}/PDC/analysis/lasso_prognosis/LASSO_enhanced_summary.csv"
lasso_sel_file = f"{base}/PDC/analysis/lasso_prognosis/LASSO_enhanced_selected_cpgs.csv"
model_cmp_file = f"{base}/PDC/analysis/model_selection_14cpg/Model_Comparison_14CpG.csv"
risk_file = f"{base}/PDC/analysis/model_selection_14cpg/BestModel_RiskScores.csv"
figure7_sum_file = f"{base}/Final_Submission_Figures/Main_Figures_14CpG/Figure7_M12_Summary.csv"
full_matrix_file = f"{base}/PDC/analysis/lasso_prognosis/LASSO_Risk_Scores.csv"

sup_digits = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")

def fmt_p(p):
    try:
        v = float(p)
    except Exception:
        return "NA"
    if not np.isfinite(v):
        return "NA"
    if v < 0.001:
        if v == 0:
            return "<0.001"
        e = int(np.floor(np.log10(v)))
        m = v / (10 ** e)
        return f"{m:.2f} × 10{str(e).translate(sup_digits)}"
    return f"{v:.3f}"

def fmt_num(v, d=3):
    try:
        x = float(v)
    except Exception:
        return "NA"
    if not np.isfinite(x):
        return "NA"
    return f"{x:.{d}f}"

def parse_p_any(v):
    if v is None:
        return np.nan
    s = str(v).strip().replace(",", "")
    if s == "" or s.lower() in {"na", "n/a", "nan"}:
        return np.nan
    m = None
    if "x 10^{" in s:
        m = __import__("re").search(r"([-+]?\d*\.?\d+)\s*x\s*10\^\{([-+]?\d+)\}", s, flags=__import__("re").I)
    if m:
        try:
            return float(m.group(1)) * (10 ** int(m.group(2)))
        except Exception:
            return np.nan
    if "×" in s and "10" in s:
        m2 = __import__("re").search(r"([-+]?\d*\.?\d+)\s*×\s*10([⁻⁰¹²³⁴⁵⁶⁷⁸⁹-]+)", s)
        if m2:
            trans = str.maketrans("⁻⁰¹²³⁴⁵⁶⁷⁸⁹", "-0123456789")
            try:
                return float(m2.group(1)) * (10 ** int(m2.group(2).translate(trans)))
            except Exception:
                return np.nan
    try:
        return float(s)
    except Exception:
        mm = __import__("re").search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", s)
        return float(mm.group(0)) if mm else np.nan

def normalize_ci_text(s):
    t = str(s).strip()
    t = t.replace("--", "–").replace(" - ", "–")
    t = re.sub(r"(?<=\d)-(?=\d)", "–", t)
    return t

def read_clean_supp_csv(path):
    rows = []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            t = line.strip().replace("\ufeff", "")
            if not t:
                continue
            low = t.lower()
            if t.startswith("\\") or "addcontentsline" in low or "vspace" in low or "scriptsize" in low or "centering" in low:
                continue
            if low.startswith("this table") or low.startswith('"this table'):
                continue
            rows.append(t)
    if len(rows) == 0:
        return pd.DataFrame()
    from io import StringIO
    text = "\n".join(rows)
    d = pd.read_csv(StringIO(text))
    # drop duplicated header rows
    if len(d) > 1:
        mask_dup_header = d.apply(lambda r: all(str(r[c]).strip() == str(c).strip() for c in d.columns), axis=1)
        d = d[~mask_dup_header].copy()
    return d.reset_index(drop=True)

def normalize_missing_df(df):
    if df is None or len(df) == 0:
        return df
    out = df.copy()
    out = out.replace(r"^\s*$", np.nan, regex=True)
    return out.fillna("NA")

# ---------- S1-S3: preserved MR module ----------
s1 = read_clean_supp_csv(s1_src_file)
s2 = read_clean_supp_csv(s2_src_file)
s3 = read_clean_supp_csv(s3_src_file)
consortium_col = next((c for c in s1.columns if str(c).strip().lower() == "consortium"), None) if len(s1.columns) > 0 else None
if consortium_col is not None:
    s1[consortium_col] = s1[consortium_col].replace(r"^\s*$", np.nan, regex=True).fillna("NA")
s2 = s2.rename(columns={"P-value": "P-value", "Heterogeneity P": "Heterogeneity P-value", "Pleiotropy P": "Pleiotropy P-value"})
for c in ["P-value", "Heterogeneity P-value", "Pleiotropy P-value"]:
    if c in s2.columns:
        s2[c] = s2[c].apply(lambda x: fmt_p(parse_p_any(x)))
if "OR (95% CI)" in s3.columns:
    s3["OR (95% CI)"] = s3["OR (95% CI)"].apply(normalize_ci_text)
if "P-value" in s3.columns:
    s3["P-value"] = s3["P-value"].apply(lambda x: fmt_p(parse_p_any(x)))
s1 = normalize_missing_df(s1)
s2 = normalize_missing_df(s2)
s3 = normalize_missing_df(s3)
s1.to_csv(f"{out_dir}/Supplementary_Table_S1.csv", index=False)
s2.to_csv(f"{out_dir}/Supplementary_Table_S2.csv", index=False)
s3.to_csv(f"{out_dir}/Supplementary_Table_S3.csv", index=False)

# ---------- S4: full colocalization ----------
coloc = pd.read_csv(coloc_file)
s4_cols = [
    "Gene_Symbol", "ENSEMBL_ID", "Chromosome", "Position_Start", "Position_End",
    "Colocalization_Score", "SMR_P_Value", "Hweis_P_Value",
    "N_Mendelian_Randomization_Studies", "Primary_Exposure_Trait", "Primary_Outcome_Trait",
]
s4_cols = [c for c in s4_cols if c in coloc.columns]
s4 = coloc[s4_cols].copy()
s4 = s4.sort_values(by=[c for c in ["Colocalization_Score", "SMR_P_Value"] if c in s4.columns], ascending=[False, True][:len([c for c in ["Colocalization_Score", "SMR_P_Value"] if c in s4.columns])])
if "Colocalization_Score" in s4.columns:
    s4["Colocalization_Score"] = s4["Colocalization_Score"].apply(lambda x: fmt_num(x, 3))
if "SMR_P_Value" in s4.columns:
    s4["SMR_P_Value"] = s4["SMR_P_Value"].apply(fmt_p)
if "Hweis_P_Value" in s4.columns:
    s4["Hweis_P_Value"] = s4["Hweis_P_Value"].apply(fmt_p)
if "N_Mendelian_Randomization_Studies" in s4.columns:
    s4["N_Mendelian_Randomization_Studies"] = pd.to_numeric(
        s4["N_Mendelian_Randomization_Studies"], errors="coerce"
    ).fillna(0).astype(int).astype(str)
s4.to_csv(f"{out_dir}/Supplementary_Table_S4.csv", index=False)

# ---------- S5: CpG annotation ----------
coef = pd.read_csv(coef_file)
bridge = pd.read_csv(bridge_file)
coef["Feature"] = coef["Feature"].astype(str)
s5 = coef.merge(bridge, left_on="Feature", right_on="CpG", how="left")
g = s5.groupby("Feature", as_index=False)["Gene"].apply(lambda x: "; ".join(sorted(set([str(v) for v in x.dropna() if str(v) != "NA"]))))
if isinstance(g, pd.DataFrame):
    g = g.rename(columns={"Gene": "Gene_Annotation"})
else:
    g = g.reset_index().rename(columns={0: "Gene_Annotation"})
s5 = coef.merge(g, on="Feature", how="left")
s5["Gene_Annotation"] = s5["Gene_Annotation"].replace("", np.nan).fillna("intergenic")
s5["Direction"] = np.where(s5["Coef"] > 0, "Risk-increasing", "Protective")
s5 = s5.rename(columns={"Feature": "CpG_ID", "Coef": "Coefficient_beta", "P": "P_value"})
s5 = s5[["CpG_ID", "Gene_Annotation", "Coefficient_beta", "HR", "P_value", "Direction"]]
s5["Coefficient_beta"] = s5["Coefficient_beta"].apply(lambda x: fmt_num(x, 3))
s5["HR"] = s5["HR"].apply(lambda x: fmt_num(x, 2))
s5["P_value"] = s5["P_value"].apply(fmt_p)
s5.to_csv(f"{out_dir}/Supplementary_Table_S5.csv", index=False)

# ---------- S6: Full M12 coefficients with CI + selection frequency ----------
mat = pd.read_csv(full_matrix_file)
keep_ids = set(pd.read_csv(risk_file)["Sample_ID"].astype(str))
mat = mat[mat["Sample_ID"].astype(str).isin(keep_ids)].copy()
features = s5["CpG_ID"].tolist()
for f in features:
    mat[f] = pd.to_numeric(mat[f], errors="coerce")
    med = mat[f].median()
    mat[f] = mat[f].fillna(med)
    sd = mat[f].std(ddof=0)
    if not np.isfinite(sd) or sd == 0:
        mat[f + "_z"] = 0.0
    else:
        mat[f + "_z"] = (mat[f] - mat[f].mean()) / sd
cox_df = pd.DataFrame({
    "time": pd.to_numeric(mat["OS_Time"], errors="coerce") / 365.25,
    "event": pd.to_numeric(mat["OS_Event"], errors="coerce"),
})
for f in features:
    cox_df[f] = mat[f + "_z"].values
cox_df = cox_df.replace([np.inf, -np.inf], np.nan).dropna()
cox_df = cox_df[(cox_df["time"] > 0) & (cox_df["event"].isin([0, 1]))].copy()
cox_df["event"] = cox_df["event"].astype(int)

cf = CoxPHFitter()
cf.fit(cox_df, duration_col="time", event_col="event")
coef_tab = cf.summary.reset_index().rename(columns={"covariate": "CpG_ID"})
coef_tab["CpG_ID"] = coef_tab["CpG_ID"].astype(str)
coef_tab["HR"] = np.exp(coef_tab["coef"])
coef_tab["HR_CI_Lower"] = np.exp(coef_tab["coef lower 95%"])
coef_tab["HR_CI_Upper"] = np.exp(coef_tab["coef upper 95%"])
coef_tab = coef_tab[["CpG_ID", "coef", "HR", "HR_CI_Lower", "HR_CI_Upper", "p"]].rename(columns={"coef": "Coefficient_beta", "p": "P_value"})

sel = pd.read_csv(lasso_sel_file)
n_models = sel["model"].nunique()
freq = sel.groupby("cpg")["model"].nunique().reset_index().rename(columns={"cpg": "CpG_ID", "model": "Selection_Count"})
freq["Selection_Frequency"] = freq["Selection_Count"] / max(n_models, 1)

s6 = s5.merge(coef_tab, on="CpG_ID", how="left", suffixes=("", "_mv"))
s6 = s6.merge(freq[["CpG_ID", "Selection_Count", "Selection_Frequency"]], on="CpG_ID", how="left")
s6["Selection_Count"] = s6["Selection_Count"].fillna(0).astype(int)
s6["Selection_Frequency"] = s6["Selection_Frequency"].fillna(0.0)
s6 = s6[[
    "CpG_ID", "Gene_Annotation", "Coefficient_beta_mv", "HR_mv", "HR_CI_Lower", "HR_CI_Upper",
    "P_value_mv", "Direction", "Selection_Count", "Selection_Frequency"
]].rename(columns={
    "Coefficient_beta_mv": "Coefficient_beta",
    "HR_mv": "HR",
    "P_value_mv": "P_value"
})
s6["Direction"] = s6["Direction"].replace({"Risk-increasing": "Risk-increasing", "Protective": "Protective"})
s6_fmt = pd.DataFrame({
    "CpG ID": s6["CpG_ID"],
    "Gene": s6["Gene_Annotation"],
    "β": s6["Coefficient_beta"].apply(lambda x: fmt_num(x, 3)),
    "HR (95% CI)": s6.apply(lambda r: f"{fmt_num(r['HR'],2)} ({fmt_num(r['HR_CI_Lower'],2)}–{fmt_num(r['HR_CI_Upper'],2)})", axis=1),
    "P-value": s6["P_value"].apply(fmt_p),
    "Direction": s6["Direction"].replace({"Risk-increasing": "Risk-increasing", "Protective": "Protective"}),
    "Selection frequency": s6["Selection_Frequency"].apply(lambda x: fmt_num(x, 3)),
})
s6_fmt.to_csv(f"{out_dir}/Supplementary_Table_S6.csv", index=False)

# ---------- S7: feature selection process ----------
sumdf = pd.read_csv(lasso_sum_file).iloc[0]
cmp = pd.read_csv(model_cmp_file).sort_values("Selection_Score").reset_index(drop=True)
top5 = cmp.head(5)[["Model", "N_CpG", "Selection_Score", "C_Index_Oriented", "KM_LogRank_P", "AUC_Mean", "AIC"]]
top5 = top5.rename(columns={"Model": "Candidate model", "N_CpG": "N CpGs", "Selection_Score": "Selection score", "C_Index_Oriented": "C-index", "KM_LogRank_P": "KM P-value", "AUC_Mean": "Mean AUC", "AIC": "AIC"})
for c in ["C-index", "Mean AUC"]:
    top5[c] = top5[c].apply(lambda x: fmt_num(x, 3))
top5["KM P-value"] = top5["KM P-value"].apply(fmt_p)
top5["AIC"] = top5["AIC"].apply(lambda x: fmt_num(x, 1))
steps = pd.DataFrame([
    ["Raw features in matrix", int(sumdf["n_cpg_after_missing_filter"]), "After missing-value filter (<=20%)", "LASSO_enhanced_summary"],
    ["After univariate Cox prefilter", int(sumdf["n_cpg_after_univariate_prefilter"]), f"P-threshold={sumdf['prefilter_p_threshold']}", "LASSO_enhanced_summary"],
    ["LASSO alpha=1 selected (lambda.min)", int(sumdf["lasso_n_selected_lambda_min"]), f"lambda.min={sumdf['lasso_lambda_min']:.6f}", "LASSO_enhanced_summary"],
    ["Elastic Net alpha=0.5 selected (lambda.min)", int(sumdf["enet_n_selected_lambda_min"]), f"lambda.min={sumdf['enet_lambda_min']:.6f}", "LASSO_enhanced_summary"],
    ["Final selected model", 12, "M12_Top12 by composite Selection_Score", "Model_Comparison_14CpG"],
], columns=["Step", "N_Features", "Details", "Source"])
s7 = {
    "process_steps": steps,
    "top_models": top5
}
steps.to_csv(f"{out_dir}/Supplementary_Table_S7_Process.csv", index=False)
top5.to_csv(f"{out_dir}/Supplementary_Table_S7_TopModels.csv", index=False)

# ---------- S8: final model performance ----------
fig7 = pd.read_csv(figure7_sum_file)

tcga = pd.read_csv(risk_file).copy()
tcga["time_year"] = pd.to_numeric(tcga["OS_Time"], errors="coerce") / 365.25
tcga["event"] = pd.to_numeric(tcga["OS_Event"], errors="coerce").astype(int)
tcga = tcga.replace([np.inf, -np.inf], np.nan).dropna(subset=["Risk_Score", "time_year", "event"])
tcga = tcga[(tcga["time_year"] > 0) & (tcga["event"].isin([0, 1]))].copy()
tcga["risk01"] = (tcga["Risk_Score"] - tcga["Risk_Score"].min()) / (tcga["Risk_Score"].max() - tcga["Risk_Score"].min() + 1e-12)
tcga["y1"] = ((tcga["event"] == 1) & (tcga["time_year"] <= 1)).astype(int)
tcga["y3"] = ((tcga["event"] == 1) & (tcga["time_year"] <= 3)).astype(int)
tcga["y5"] = ((tcga["event"] == 1) & (tcga["time_year"] <= 5)).astype(int)

cidx = concordance_index(tcga["time_year"].values, -tcga["Risk_Score"].values, tcga["event"].values)
cidx = max(cidx, 1 - cidx)

cox_sd = pd.DataFrame({
    "time": tcga["time_year"].values,
    "event": tcga["event"].values,
    "risk_sd": (tcga["Risk_Score"].values - tcga["Risk_Score"].mean()) / (tcga["Risk_Score"].std(ddof=0) + 1e-12),
})
cf_sd = CoxPHFitter()
cf_sd.fit(cox_sd, duration_col="time", event_col="event")
hr_sd = np.exp(cf_sd.params_["risk_sd"])
l_sd = np.exp(cf_sd.confidence_intervals_.loc["risk_sd", "95% lower-bound"])
u_sd = np.exp(cf_sd.confidence_intervals_.loc["risk_sd", "95% upper-bound"])
p_sd = cf_sd.summary.loc["risk_sd", "p"]

def calibration_stats(d, ycol, bins=8):
    tmp = d[["risk01", ycol]].copy()
    tmp["bin"] = pd.qcut(tmp["risk01"], q=bins, labels=False, duplicates="drop")
    gp = tmp.groupby("bin").agg(pred=("risk01", "mean"), obs=(ycol, "mean")).reset_index()
    if len(gp) < 2:
        return np.nan, np.nan
    m = np.polyfit(gp["pred"], gp["obs"], 1)
    slope = float(m[0])
    r2 = float(np.corrcoef(gp["pred"], gp["obs"])[0, 1] ** 2)
    return slope, r2

slope3, r2_3 = calibration_stats(tcga, "y3")
slope5, r2_5 = calibration_stats(tcga, "y5")

auc1 = roc_auc_score(tcga["y1"], tcga["risk01"])
auc3 = roc_auc_score(tcga["y3"], tcga["risk01"])
auc5 = roc_auc_score(tcga["y5"], tcga["risk01"])

ths = np.arange(0.01, 0.51, 0.01)
y = tcga["y3"].values
p = tcga["risk01"].values
n = len(y)
nb_model = []
nb_all = []
for t in ths:
    pred = p >= t
    tp = np.sum((pred == 1) & (y == 1))
    fp = np.sum((pred == 1) & (y == 0))
    nb_model.append((tp / n) - (fp / n) * (t / (1 - t)))
    ev = y.mean()
    nb_all.append(ev - (1 - ev) * (t / (1 - t)))
nb_model = np.array(nb_model)
nb_all = np.array(nb_all)
dca_pos = ths[nb_model > np.maximum(nb_all, 0)]
dca_range = f"{dca_pos.min():.2f}–{dca_pos.max():.2f}" if len(dca_pos) > 0 else "NA"

perf_rows = []
for _, r in fig7.iterrows():
    cohort = r["Cohort"]
    endpoint = "OS" if "OS" in cohort or "TCGA" in cohort else "RFS"
    row = {
        "Dataset": "TCGA (training)" if cohort == "TCGA_train" else ("GSE39279 (external)" if cohort == "GSE39279_RFS" else "GSE30219 (external)"),
        "Endpoint": "Overall survival (OS)" if endpoint == "OS" else "Recurrence-free survival (RFS)",
        "N": int(r["N"]),
        "KM_P_value": r["LogRank_P"],
        "Cutoff_Method": "Median (training cohort)" if str(r["Cutoff_Method"]) == "median_from_training" else "Maxstat-derived optimal cutoff",
        "Cutoff_Value": r["Cutoff_Value"],
        "C_index": np.nan,
        "AUC_3Y": np.nan,
        "Calibration_Slope_3Y": np.nan,
        "Calibration_R2_3Y": np.nan,
        "HR_per_SD": np.nan,
        "HR_per_SD_P": np.nan,
    }
    if cohort == "TCGA_train":
        row["C_index"] = cidx
        row["AUC_3Y"] = auc3
        row["Calibration_Slope_3Y"] = slope3
        row["Calibration_R2_3Y"] = r2_3
        row["HR_per_SD"] = f"{hr_sd:.2f} ({l_sd:.2f}–{u_sd:.2f})"
        row["HR_per_SD_P"] = p_sd
    perf_rows.append(row)
s8 = pd.DataFrame(perf_rows)
s8_fmt = s8.copy()
s8_fmt["N"] = s8_fmt["N"].astype(int).astype(str)
s8_fmt["KM P-value"] = s8_fmt["KM_P_value"].apply(fmt_p)
s8_fmt["Cutoff value"] = s8_fmt["Cutoff_Value"].apply(lambda x: fmt_num(x, 3))
s8_fmt["C-index"] = s8_fmt["C_index"].apply(lambda x: fmt_num(x, 3))
s8_fmt["AUC (3-year)"] = s8_fmt["AUC_3Y"].apply(lambda x: fmt_num(x, 3))
s8_fmt["Calibration slope (3Y)"] = s8_fmt["Calibration_Slope_3Y"].apply(lambda x: fmt_num(x, 3))
s8_fmt["Calibration R² (3Y)"] = s8_fmt["Calibration_R2_3Y"].apply(lambda x: fmt_num(x, 3))
s8_fmt["HR (per SD, 95% CI)"] = s8_fmt["HR_per_SD"].fillna("NA")
s8_fmt["P-value (HR per SD)"] = s8_fmt["HR_per_SD_P"].apply(fmt_p)
s8_fmt = s8_fmt.rename(columns={"Cutoff_Method": "Cutoff method"})
s8_fmt = s8_fmt[[
    "Dataset", "Endpoint", "N", "HR (per SD, 95% CI)", "P-value (HR per SD)",
    "KM P-value", "C-index", "AUC (3-year)", "Calibration slope (3Y)", "Calibration R² (3Y)",
    "Cutoff method", "Cutoff value"
]]
s8_fmt.to_csv(f"{out_dir}/Supplementary_Table_S8.csv", index=False)

# ---------- S9: diagnostics summary ----------
s9 = pd.DataFrame([
    ["DCA threshold range (3-year)", dca_range, "Decision curve analysis", "TCGA M12"],
    ["AUC (1-year)", auc1, "Time-dependent ROC", "TCGA M12"],
    ["AUC (3-year)", auc3, "Time-dependent ROC", "TCGA M12"],
    ["AUC (5-year)", auc5, "Time-dependent ROC", "TCGA M12"],
    ["Calibration slope (3-year)", slope3, "Calibration", "TCGA M12"],
    ["Calibration R² (3-year)", r2_3, "Calibration", "TCGA M12"],
    ["Calibration slope (5-year)", slope5, "Calibration", "TCGA M12"],
    ["Calibration R² (5-year)", r2_5, "Calibration", "TCGA M12"],
    ["C-index", cidx, "Concordance", "TCGA M12"],
], columns=["Metric", "Value", "Method", "Cohort"])
s9["Value"] = s9["Value"].apply(lambda x: fmt_num(x, 3) if isinstance(x, (int, float, np.floating)) else x)
s9.to_csv(f"{out_dir}/Supplementary_Table_S9.csv", index=False)

# ---------- Rebuilt supplementary markdown ----------
md = f"""\n**Supplementary Tables (Rebuilt for M12-consistent submission)**\n\n**Module 1: Causal inference (preserved)**\n- Table S1A: Core MR results (unchanged)\n- Table S1B: MR sensitivity analyses (unchanged)\n- Table S2: MVMR results (unchanged)\n- Table S3: Extended MVMR analyses (unchanged)\n\n**Module 2: Gene mapping / mechanism expansion**\n- Table S4: Full colocalization results\n- Table S5: CpG-gene annotation for the M12 model\n\n**Module 3: Model construction (M12)**\n- Table S6: Full M12 model coefficients (β, HR, 95%CI, P, selection frequency)\n- Table S7: Feature-selection process summary (steps + top candidate models)\n\n**Module 4: Model validation (M12)**\n- Table S8: Final M12 model performance summary across cohorts\n- Table S9: Diagnostic metrics (DCA + ROC + calibration)\n\n---\n\n**Generated files in** `{out_dir}`:\n- Supplementary_Table_S4.csv\n- Supplementary_Table_S5.csv\n- Supplementary_Table_S6.csv\n- Supplementary_Table_S7_Process.csv\n- Supplementary_Table_S7_TopModels.csv\n- Supplementary_Table_S8.csv\n- Supplementary_Table_S9.csv\n\n**Primary model anchor**: M12_Top12\n"""

with open(f"{out_dir}/Supplementary_Tables_Final_updated.md", "w", encoding="utf-8") as f:
    f.write(md)

print("Rebuilt supplementary tables S4-S9 with M12 consistency.")
